!
!***********************************************************************
!
      module Coupling_transform_cfg_LSLK
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: main_cfg_lslk
      public  :: count_nr_of_csfs_LK
      public  :: delete_cfg_expansions
      private :: form_list_nomach_csfs
      private :: form_csfs_LK
      private :: matrix_LS_LK
      private :: dotransform_LSLK
      private :: print_cfg_LSLK
      private :: equivalent_LSLK
!-----------------------------------------------
!   D e f i n i t i o n s
!-----------------------------------------------
      type::Ji_lists
         integer::nr_of_csf !serial number of csf_LS in csfs_LS
         type(list),dimension(:),pointer::Ji !i=1..nr_of_subc
      end type Ji_lists
!
      type:: cfg_Ji_lists
         integer::nr_of_subc
         integer::nr_of_nomach_csfs
         type(Ji_lists),dimension(:),pointer::csfs_Ji !i=1..nr_of_nomach_csfs
      end type cfg_Ji_lists
!
      type(expansion), public :: expansion_cfg_LS
      type(expansion), public :: expansion_cfg_LK
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_cfg_lslk(print_level)
!                                                                      *
!     This is managerial subroutine                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in) :: print_level
      integer::error
      integer::itype
      expansion_cfg_LK%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=0
      call form_csfs_LK(itype)
      call dotransform_LSLK
      if(print_level.gt.1) call print_cfg_LSLK(2) ! 2-means print LS and LK expansions
      deallocate(nomach_csfs_ls%items, STAT = error)
      end subroutine main_cfg_lslk
!
!***********************************************************************
!                                                                      *
      subroutine count_nr_of_csfs_LK(nr_of_csfs)
!                                                                      *
!     This subroutine counts the number of oneconfigurational          *
!     expansion in LK coupling                                         *
!     (corresponding to the one in the JJ coupling "expansion_cfg_LS") *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(out):: nr_of_csfs
      integer::error
      integer::itype
      if(expansion_cfg_LS%csfs(1)%nosubc.lt.2) then
         write(*,*)                                                    &
            'ERROR at subroutine count_nr_of_csfs_LK: nosubc.lt.2'
         write(*,*)                                                    &
            '(you can not use LK coupling with less than 2 subshells)'
         write(*,*)'program is terminated'
         stop
      end if
      expansion_cfg_LK%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=1
      call form_csfs_LK(itype)
      nr_of_csfs = expansion_cfg_LK%size
      if(associated(nomach_csfs_ls%items))                             &
                             deallocate(nomach_csfs_ls%items,STAT=error)
      end subroutine count_nr_of_csfs_LK
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_expansions
!                                                                      *
!     This subroutine deallocates the arrays of                        *
!     oneconfigurational expansions "expansion_cfg_LK" and             *
!     "expansion_cfg_LS"                                               *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      integer::icsf
!
!write(iwrite_log,'(6x,a21)') 'delete_cfg_expansions'
!
!delete expansion_cfg_LK
      if(associated(expansion_cfg_LK%coeffs))                          &
                                      deallocate(expansion_cfg_LK%csfs)
      if(associated(expansion_cfg_LK%csfs)) then
         do icsf=1, expansion_cfg_LK%size
            if(associated(expansion_cfg_LK%csfs(icsf)%subc_cfg))       &
                        deallocate(expansion_cfg_LK%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_LK%csfs(icsf)%subc))           &
                            deallocate(expansion_cfg_LK%csfs(icsf)%subc)
            if(associated(expansion_cfg_LK%csfs(icsf)%iM1))            &
                             deallocate(expansion_cfg_LK%csfs(icsf)%iM1)
            if(associated(expansion_cfg_LK%csfs(icsf)%iM2))            &
                             deallocate(expansion_cfg_LK%csfs(icsf)%iM2)
            if(associated(expansion_cfg_LK%csfs(icsf)%iJ))            &
                             deallocate(expansion_cfg_LK%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_LK%csfs)
      end if
      expansion_cfg_LK%size=0
      expansion_cfg_LK%nr_of_state=0
!
!delete expansion_cfg_LS
      if(associated(expansion_cfg_LS%coeffs))                          &
                                       deallocate(expansion_cfg_LS%csfs)
      if(associated(expansion_cfg_LS%csfs)) then
         do icsf=1, expansion_cfg_LS%size
            if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg))       &
                        deallocate(expansion_cfg_LS%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_LS%csfs(icsf)%subc))           &
                            deallocate(expansion_cfg_LS%csfs(icsf)%subc)
            if(associated(expansion_cfg_LS%csfs(icsf)%iM1))            &
                             deallocate(expansion_cfg_LS%csfs(icsf)%iM1)
            if(associated(expansion_cfg_LS%csfs(icsf)%iM2))            &
                             deallocate(expansion_cfg_LS%csfs(icsf)%iM2)
            if(associated(expansion_cfg_LS%csfs(icsf)%iJ))            &
                             deallocate(expansion_cfg_LS%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_LS%csfs)
      end if
      expansion_cfg_LS%size=0
      expansion_cfg_LS%nr_of_state=0
!
!write(iwrite_log,'(6x,a25)') 'end delete_cfg_expansions'
!
      end subroutine delete_cfg_expansions
!
!***********************************************************************
!                                                                      *
      function equivalent_LSLK(csf1, csf2)   result(rez)
!                                                                      *
!     This function defines the "equivalency" of                       *
!     two csfs in LS coupling for the formation of LK csfs             *
!     (for notion of the "equivalency" see the                         *
!     description in the program)                                      *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(csf_LS) :: csf1, csf2
      logical :: rez
      integer::isubc
      rez = .true.
      if(csf1%nosubc.lt.2.or.csf2%nosubc.lt.2) then
         write(*,*)'ERROR at subroutine equivalent_LSLK: ',            &
                   'csf1%nosubc.lt.2.or.csf2%nosubc.lt.2'
         write(*,*)'(you can not use LK coupling with less than 2 ',   &
                   'subshells)'
         write(*,*)'program will be terminated'
         stop
      end if
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc-1
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%iM1(isubc)==csf2%iM1(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%iM2(isubc)==csf2%iM2(isubc)) then
               rez=.false.
               exit
            end if
         end do
         isubc=csf1%nosubc
         if(rez) then
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc))       &
                                                           rez=.false.
            if(.not.csf1%subc(isubc)==csf2%subc(isubc)) rez=.false.
            if(.not.csf1%iM1(isubc)==csf2%iM1(isubc)) rez=.false.
         end if
      else
         rez=.false.
      end if
      end function equivalent_LSLK
!
!***********************************************************************
!                                                                      *
      subroutine form_list_nomach_csfs
!                                                                      *
!     This subroutine form the list of serial numbers                  *
!     in "expansion_cfg_LS" of "nonequivalent" LS coupling csfs        *
!     (for notion of the "equivalency" see the                         *
!     description in the program)                                      *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer:: inomach_counter
      integer, dimension(:), pointer :: temp_list
      integer::icsf_LS, inew
      logical:: new_csf
!
      inomach_counter=0
      allocate(temp_list(expansion_cfg_LS%size))
!
      do icsf_LS=1,expansion_cfg_LS%size,1
         new_csf = .true.
         do inew = 1, inomach_counter, 1
            if(equivalent_LSLK(expansion_cfg_LS%csfs(icsf_LS),         &
              expansion_cfg_LS%csfs(temp_list(inew)))) new_csf = .false.
         end do
         if(new_csf) then
            inomach_counter = inomach_counter + 1
            if(inomach_counter.gt.expansion_cfg_LS%size) then
               write(*,*) 'stop at subroutine define_nomach: ',        &
                  'inomach_counter.gt.expansion_cfg_LS%size'
               stop
            end if
            temp_list(inomach_counter)= icsf_LS
         end if
      end do
!
      nomach_csfs_LS%list_size = inomach_counter
!
      allocate(nomach_csfs_LS%items(nomach_csfs_LS%list_size))
!
      do icsf_LS=1, nomach_csfs_LS%list_size, 1
         nomach_csfs_LS%items(icsf_LS) = temp_list(icsf_LS)
      end do
!
      deallocate(temp_list)
!----- print down data -----------
!write(iwrite_log, *)'Nonmaching csfs_LS:'
!write(iwrite_log, *)'       Nr.   icsf_nr '
!do icsf_LS=1, nomach_csfs_LS%list_size, 1
!	 write(iwrite_log, '(4x,i3,2x,i3)') icsf_LS, nomach_csfs_LS%items(icsf_LS)
!end do !icsf_LS
!write(iwrite_log, *)'  '
!
      end subroutine form_list_nomach_csfs
!
!***********************************************************************
!                                                                      *
      subroutine form_csfs_LK(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in LK coupling "expansion_cfg_LK",                     *
!     correponding to the one in LS coupling "expansion_cfg_LS"        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
!type(cfg_Ji_lists),intent(in)::cfg_Ji_structure
      integer:: itype !defines the type of execution
                      !itype=1 - to find the number
                      !of csfs_LK only,
                      !itype.ne.1 - perform full
                      !calculation
      integer :: inr_of_csfs_LK   !, icsf_nr
      integer :: icsf_LS, icsf_LK
      integer :: K, K_min, K_max, J, J_min,J_max
      integer :: isubc_aviable
      integer :: inosubc, isubc
      integer :: iJ_total
      integer :: inr_of_csf_in_expansion_LS
!
      isubc_aviable = 2
      if(expansion_cfg_LS%csfs(1)%nosubc.gt.isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_LK ',    &
            'module transform_lsjj: cfg_Ji_structure%nr_of_subc.gt.',  &
            'isubc_aviable'
         stop
      end if
!
!write(*,*) '      subroutine form_csfs_LK'
!write(iwrite_log,*) '      subroutine form_csfs_LK'
!
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
!
      call define_number_of_csfs_LK(inr_of_csfs_LK)
!
      expansion_cfg_LK%size=inr_of_csfs_LK
      if(itype.ne.1) then
         allocate(expansion_cfg_LK%csfs(expansion_cfg_LK%size))
         allocate(expansion_cfg_LK%coeffs(expansion_cfg_LK%size))
         icsf_LK=0
         do icsf_LS=1, nomach_csfs_ls%list_size,1
            inr_of_csf_in_expansion_LS = nomach_csfs_ls%items(icsf_LS)
            inosubc = expansion_cfg_LS%csfs(                           &
                                     inr_of_csf_in_expansion_LS)%nosubc
            K_min = abs(expansion_cfg_LS%csfs(                         &
               inr_of_csf_in_expansion_LS)%iM1(inosubc) -              &
               expansion_cfg_LS%csfs(inr_of_csf_in_expansion_LS)%      &
               iM2(inosubc-1))
            K_max=abs(expansion_cfg_LS%csfs(                           &
               inr_of_csf_in_expansion_LS)%iM1(inosubc) +              &
               expansion_cfg_LS%csfs(inr_of_csf_in_expansion_LS)%      &
               iM2(inosubc-1))
            do K= K_min, K_max, 2
               J_min = abs(K - expansion_cfg_LS%csfs(                  &
                  inr_of_csf_in_expansion_LS)%subc(inosubc)%iS)
               J_max=abs(K + expansion_cfg_LS%csfs(                    &
                  inr_of_csf_in_expansion_LS)%subc(inosubc)%iS)
               do J = J_min, J_max, 2
                  if(J.eq.iJ_total) then
                     icsf_LK = icsf_LK + 1
                     if(icsf_LK.gt.expansion_cfg_LK%size) then
                        stop
                        write(8,*) 'stop at subroutine form_csfs_LK: ',&
                           'icsf_nr.gt.expansion_cfg_LK%size'
                     end if
                     expansion_cfg_LK%csfs(icsf_LK)%nosubc=inosubc
                     allocate(expansion_cfg_LK%csfs(icsf_LK)%subc_cfg( &
                        expansion_cfg_LK%csfs(icsf_LK)%nosubc))
                     allocate(expansion_cfg_LK%csfs(icsf_LK)%subc(     &
                        expansion_cfg_LK%csfs(icsf_LK)%nosubc))
                     allocate(expansion_cfg_LK%csfs(icsf_LK)%iM1(      &
                        expansion_cfg_LK%csfs(icsf_LK)%nosubc))
                     allocate(expansion_cfg_LK%csfs(icsf_LK)%iM2(      &
                        expansion_cfg_LK%csfs(icsf_LK)%nosubc))
                     allocate(expansion_cfg_LK%csfs(icsf_LK)%iJ(       &
                        expansion_cfg_LK%csfs(icsf_LK)%nosubc))
                     expansion_cfg_LK%csfs(icsf_LK)%iJ(                &
                     expansion_cfg_LK%csfs(icsf_LK)%nosubc) = 0
                     do isubc=1,expansion_cfg_LK%csfs(icsf_LK)%        &
                                                            nosubc-1,1
                        expansion_cfg_LK%csfs(icsf_LK)%iJ(isubc) = 0
                        expansion_cfg_LK%csfs(icsf_LK)%subc_cfg(isubc) &
                           = expansion_cfg_LS%csfs(                    &
                           inr_of_csf_in_expansion_LS)%subc_cfg(isubc)
                        expansion_cfg_LK%csfs(icsf_LK)%subc(isubc) =   &
                           expansion_cfg_LS%csfs(                      &
                           inr_of_csf_in_expansion_LS)%subc(isubc)
                        expansion_cfg_LK%csfs(icsf_LK)%iM1(isubc) =    &
                           expansion_cfg_LS%csfs(                      &
                           inr_of_csf_in_expansion_LS)%iM1(isubc)
                        expansion_cfg_LK%csfs(icsf_LK)%iM2(isubc) =    &
                           expansion_cfg_LS%csfs(                      &
                           inr_of_csf_in_expansion_LS)%iM2(isubc)
                     end do
                     isubc = expansion_cfg_LK%csfs(icsf_LK)%nosubc
                     expansion_cfg_LK%csfs(icsf_LK)%subc_cfg(isubc) =  &
                        expansion_cfg_LS%csfs(                         &
                        inr_of_csf_in_expansion_LS)%subc_cfg(isubc)
                     expansion_cfg_LK%csfs(icsf_LK)%subc(isubc) =      &
                        expansion_cfg_LS%csfs(                         &
                        inr_of_csf_in_expansion_LS)%subc(isubc)
                     expansion_cfg_LK%csfs(icsf_LK)%iM1(isubc) =       &
                        expansion_cfg_LS%csfs(                         &
                        inr_of_csf_in_expansion_LS)%iM1(isubc)
                     expansion_cfg_LK%csfs(icsf_LK)%iM2(isubc)= K
!
                  end if
               end do
            end do
         end do
!
      end if
!
!write(*,*) '      end subroutine form_csfs_LK'
!write(iwrite_log,*) '      end subroutine form_csfs_LK'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine  define_number_of_csfs_LK(irez)
!                                                                      *
!     This subroutine defines the number of csfs in LK coupling        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
         implicit none
         integer, intent(out)::irez
         integer::K, K_min, K_max, J, J_min,J_max
         integer::isubc_aviable
         integer::icsf
         integer :: icsf_nr_in_expansion_cfg_LS, inosubc
         isubc_aviable = 2
         if(expansion_cfg_LS%csfs(1)%nosubc .gt. isubc_aviable) then
           write(*,*) 'STOP at subroutine define_number_of_csfs_LK ',  &
              'module transform_lsjj: cfg_Ji_structure%nr_of_subc.gt.',&
              'isubc_aviable'
           stop
         end if
         irez=0
         do icsf=1, nomach_csfs_ls%list_size,1
            icsf_nr_in_expansion_cfg_LS = nomach_csfs_ls%items(icsf)
            inosubc = expansion_cfg_LS%csfs(                           &
               icsf_nr_in_expansion_cfg_LS)%nosubc
            K_min = abs(expansion_cfg_LS%csfs(                         &
               icsf_nr_in_expansion_cfg_LS)%iM1(inosubc) -             &
               expansion_cfg_LS%csfs(icsf_nr_in_expansion_cfg_LS)%     &
               iM2(inosubc-1))
            K_max = abs(expansion_cfg_LS%csfs(                         &
               icsf_nr_in_expansion_cfg_LS)%iM1(inosubc) +             &
               expansion_cfg_LS%csfs(icsf_nr_in_expansion_cfg_LS)%     &
               iM2(inosubc-1))
            do K= K_min, K_max, 2
               J_min = abs(K - expansion_cfg_LS%csfs(                  &
                  icsf_nr_in_expansion_cfg_LS)%subc(inosubc)%iS)
               J_max = abs(K + expansion_cfg_LS%csfs(                  &
                  icsf_nr_in_expansion_cfg_LS)%subc(inosubc)%iS)
               do J = J_min, J_max, 2
                  if(J.eq.iJ_total) irez = irez + 1
               end do
            end do
         end do
         end subroutine  define_number_of_csfs_LK
      end subroutine form_csfs_LK
!
!***********************************************************************
!                                                                      *
      subroutine matrix_LS_LK(icsf_LS,icsf_LK,rez)
!                                                                      *
!     This subroutine calculates the transformation                    *
!     matrix element between two csfs  (in LS and LK couplings)        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer,intent(in):: icsf_LS,icsf_LK
      real(kind=dp),intent(out)::rez
      integer::K, iSi, iS_i, iS_im1, iL_i, J
      integer::nosubc !number of subshells
      nosubc=expansion_cfg_LS%csfs(1)%nosubc
      rez = ZERO_dp
      if(equivalent_LSLK(expansion_cfg_LS%csfs(icsf_LS),               &
         expansion_cfg_LK%csfs(icsf_LK))) then
         iSi=expansion_cfg_LS%csfs(icsf_LS)%subc(nosubc)%iS
         iL_i=expansion_cfg_LS%csfs(icsf_LS)%iM1(nosubc)
         iS_i=expansion_cfg_LS%csfs(icsf_LS)%iM2(nosubc)
         iS_im1=expansion_cfg_LS%csfs(icsf_LS)%iM2(nosubc-1)
         K=expansion_cfg_LK%csfs(icsf_LK)%iM2(nosubc)
         J=states%states(expansion_cfg_LS%nr_of_state)%J
         call LKPER(0,0,iS_im1,iSi,iL_i,iS_i,J,iL_i,K,rez)
      end if
      end subroutine matrix_LS_LK
!
!***********************************************************************
!                                                                      *
      subroutine dotransform_LSLK
!                                                                      *
!     This subroutine calculates the weights of the                    *
!     expansions in the LK coupling                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp)::coeff_LK,coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_LK
      melement=0.
      sum_of_state = 0.
      do icsf_LK=1,expansion_cfg_LK%size
         coeff_LK = 0.
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_LK(icsf_LS,icsf_LK,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_LK = coeff_LK + coeff_LS*melement
         end do
         if((dabs(coeff_LK)-dp_coeff_precision).gt.ONE_dp) then
            write(*,*)'possible error at subroutine dotransform_LSLK:',&
               ' coeff_LK=',coeff_LK
            write(iwrite_log,*)'possible error at subroutine ',        &
               'dotransform_LSLK: coeff_LK=',coeff_LK
         end if
         expansion_cfg_LK%coeffs(icsf_LK)= coeff_LK
         sum_of_state = sum_of_state + coeff_LK*coeff_LK
      end do
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_LSLK: ',  &
            'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
            'dotransform_LSLK: sum_of_state=',sum_of_state
      end if
      end subroutine dotransform_LSLK
!
!***********************************************************************
!                                                                      *
      subroutine print_cfg_LSLK(itype)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the unit "iwrite_cfg_expansions_LK"                *
!     (see module "Coupling_constants")                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in) :: itype
      integer             :: icsf, isubc, j
      character(len=1)    :: CVAL
      character(len=4)    :: JVAL
!
!     LS - Coupling
!
      if(itype.gt.0) then
         write(iwrite_cfg_expansions_LK,*)  &
         '-------------------------------------'
         write(iwrite_cfg_expansions_LK,*)  &
         'state Nr.',expansion_cfg_LS%nr_of_state
         write(iwrite_cfg_expansions_LK,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_LS%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LS%nr_of_state)%energy
         write(iwrite_cfg_expansions_LK,'(3x,a29,i2)')               &
         'expansion size (LS coupling): ', expansion_cfg_LS%size
         write(iwrite_cfg_expansions_LK,*)''
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1,expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_LK,*)'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_LK,                   &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf,  &
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%in,       &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%iN_big,   &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_LK,                   &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
                 'To many coupled shells in Coupling_transform_cfg_LSLK'
                  end if
               else
                  write(iwrite_cfg_expansions_LK,'(9x,a51)')   &
                  'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LS%csfs(icsf)%subc)  &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM1) &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_LK,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_LS%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LS%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_LK,                 &
                     '(12x,7(1x,i2,a1,i1))')                         &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,      &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL),&
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,       &
                      j=1,2),                                        &
                      expansion_cfg_LS%csfs(icsf)%iM2(2)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(2))
                  end if
               else
                   write(iwrite_cfg_expansions_LK,'(9x,a63)') &
       'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_LK,*)' '
            end do
         else
            write(iwrite_cfg_expansions_LK,*) &
            'expansion_cfg_LS%csfs NOT associated'
         end if
      end if
!
!     LK - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_LK,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_LK,*) &
         'state Nr.',expansion_cfg_LK%nr_of_state
         write(iwrite_cfg_expansions_LK,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_LK%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LK%nr_of_state)%energy
         write(iwrite_cfg_expansions_LK,'(3x,a29,i2)')               &
         'expansion size (LK coupling): ', expansion_cfg_LK%size
         write(iwrite_cfg_expansions_LK,*)''
         if(associated(expansion_cfg_LK%csfs)) then
            do icsf=1,expansion_cfg_LK%size
               if(associated(expansion_cfg_LK%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_LK,*)'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_LK,                   &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf,  &
                     expansion_cfg_LK%csfs(icsf)%subc_cfg(1)%in,       &
                    CVAL(1,expansion_cfg_LK%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LK%csfs(icsf)%subc_cfg(1)%iN_big,   &
                     'coeff:',expansion_cfg_LK%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_LK,                   &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LK%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LK%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LK%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_LK%coeffs(icsf)
                  end if
               else
                  write(iwrite_cfg_expansions_LK,'(9x,a51)') &
                  'expansion_cfg_LK%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LK%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_LK%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_LK%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_LK,                &
                     '(12x,2(1x,i2,a1,i1))')                        &
                     expansion_cfg_LK%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LK%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LK%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_LK,                &
                     '(12x,2(1x,i2,a1,i1),3x,a1," [",a4,"]",a4)')   &
                     (expansion_cfg_LK%csfs(icsf)%subc(j)%iS+1,     &
                     CVAL(2,expansion_cfg_LK%csfs(icsf)%subc(j)%iL),&
                     expansion_cfg_LK%csfs(icsf)%subc(j)%inr,       &
                     j=1,2),                                        &
                     CVAL(2,expansion_cfg_LK%csfs(icsf)%iM1(2)),    &
                     JVAL(expansion_cfg_LK%csfs(icsf)%iM2(2)),      &
                     JVAL(states%states(expansion_cfg_LK%nr_of_state)%J)
                  end if
               else
                  write(iwrite_cfg_expansions_LK,'(9x,a63)') &
       'expansion_cfg_LK%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_LK,*)' '
            end do
         else
            write(iwrite_cfg_expansions_LK,*) &
            'expansion_cfg_LK%csfs NOT associated'
         end if
            write(iwrite_cfg_expansions_LK,*) &
            '-------------------------------------'
      end if
      end subroutine print_cfg_LSLK
!
      end module Coupling_transform_cfg_LSLK
