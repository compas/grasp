      module Coupling_transform_cfg_JJJK
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
!
      public :: main_cfg_JJJK
      public :: count_nr_of_csfs_JK
      public :: delete_cfg_expansions
      private:: form_list_nomach_csfs
      private:: form_csfs_JK
      private:: matrix_JJ_JK
      private:: dotransform_JJJK
      private:: print_cfg_JJJK
      private :: equivalent_JJJK
!
      type::Ji_lists
!integer::nr_of_subc
         integer::nr_of_csf !serial number of csf_JJ in csfs_JJ
         type(list),dimension(:),pointer::Ji !i=1..nr_of_subc
      end type Ji_lists
!
      type:: cfg_Ji_lists
         integer::nr_of_subc
         integer::nr_of_nomach_csfs
         type(Ji_lists),dimension(:),pointer::csfs_Ji !i=1..nr_of_nomach_csfs
      end type cfg_Ji_lists
!
      type(expansion), public :: expansion_cfg_JJ
      type(expansion), public :: expansion_cfg_JK
!
      type(list), private ::nomach_csfs_JJ
      type(cfg_Ji_lists), private::cfg_Ji_structure
!
      contains
!
      subroutine main_cfg_JJJK(print_level)
!subroutine main_cfg_JJjj(in_expansion_JJ, expansion_JK)
!--------------------------------------------------------------------
! This is managerial subroutine
!--------------------------------------------------------------------
      implicit none
      integer, intent(in):: print_level
      integer::error
      integer::itype
!
      if(expansion_cfg_JJ%csfs(1)%nosubc.lt.2) then
         write(*,*)'ERROR at subroutine main_cfg_JJJK: nosubc.lt.2'
         write(*,*)                                                    &
            '(you can not use JK coupling with less than 2 subshelJJ)'
         write(*,*)'program is terminated'
         stop
      end if
!
!write(*,*) '    subroutine main_cfg_JJJK'
!
      expansion_cfg_JK%nr_of_state=expansion_cfg_JJ%nr_of_state
!
      call form_list_nomach_csfs
!
      itype=0
      call form_csfs_JK(itype)
!
      call dotransform_JJJK
!
      if(print_level.gt.1) call print_cfg_JJJK(2) ! 2-means print JJ and JK expansions
!
      deallocate(nomach_csfs_JJ%items, STAT = error)
!
!call delete_cfg_expansion
!
!write(*,*) '    subroutine main_cfg_JJjj'
      end subroutine main_cfg_JJJK
!
      subroutine count_nr_of_csfs_JK(nr_of_csfs)
!--------------------------------------------------------------------
! This subroutine counts the number of oneconfigurational
! expansion in JK coupling
! (corresponding to the one in the JJ coupling "expansion_cfg_JJ")
!--------------------------------------------------------------------
      implicit none
      integer, intent(out):: nr_of_csfs
      integer::error
      integer::itype
!
!write(*,*) '    subroutine count_nr_of_csfs_JK'
!		call print_cfg_JJJK(1)
      if(expansion_cfg_JJ%csfs(1)%nosubc.lt.2) then
         write(*,*)'ERROR at count_nr_of_csfs_JK: nosubc.lt.2'
         write(*,*)                                                    &
            '(you can not use JK coupling with less than 2 subshelJJ)'
         write(*,*)'program is terminated'
         stop
      end if
!
      expansion_cfg_JK%nr_of_state=expansion_cfg_JJ%nr_of_state
!
      call form_list_nomach_csfs
!
!
      itype=1
      call form_csfs_JK(itype)
!
      nr_of_csfs = expansion_cfg_JK%size
!
      if(associated(nomach_csfs_JJ%items))                             &
                         deallocate(nomach_csfs_JJ%items, STAT = error)
!
!write(*,*) '    end subroutine count_nr_of_csfs_JK'
      end subroutine count_nr_of_csfs_JK
!
      subroutine delete_cfg_expansions
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of
! oneconfigurational expansions "expansion_cfg_JK" and
! "expansion_cfg_JJ"
!--------------------------------------------------------------------
      integer::icsf
!
!
!delete expansion_cfg_JK
      if(associated(expansion_cfg_JK%coeffs))                          &
                                      deallocate(expansion_cfg_JK%csfs)
      if(associated(expansion_cfg_JK%csfs)) then
         do icsf=1, expansion_cfg_JK%size
            if(associated(expansion_cfg_JK%csfs(icsf)%subc_cfg))       &
                       deallocate(expansion_cfg_JK%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_JK%csfs(icsf)%subc))           &
                           deallocate(expansion_cfg_JK%csfs(icsf)%subc)
            if(associated(expansion_cfg_JK%csfs(icsf)%iM1))            &
                            deallocate(expansion_cfg_JK%csfs(icsf)%iM1)
            if(associated(expansion_cfg_JK%csfs(icsf)%iM2))            &
                            deallocate(expansion_cfg_JK%csfs(icsf)%iM2)
            if(associated(expansion_cfg_JK%csfs(icsf)%iJ))            &
                            deallocate(expansion_cfg_JK%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_JK%csfs)
      end if
      expansion_cfg_JK%size=0
      expansion_cfg_JK%nr_of_state=0
!
!delete expansion_cfg_JJ
      if(associated(expansion_cfg_JJ%coeffs))                          &
                                     deallocate(expansion_cfg_JJ%csfs)
      if(associated(expansion_cfg_JJ%csfs)) then
         do icsf=1, expansion_cfg_JJ%size
            if(associated(expansion_cfg_JJ%csfs(icsf)%subc_cfg))       &
                        deallocate(expansion_cfg_JJ%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_JJ%csfs(icsf)%subc))           &
                            deallocate(expansion_cfg_JJ%csfs(icsf)%subc)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iM1))            &
                             deallocate(expansion_cfg_JJ%csfs(icsf)%iM1)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iM2))            &
                             deallocate(expansion_cfg_JJ%csfs(icsf)%iM2)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iJ))            &
                             deallocate(expansion_cfg_JJ%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_JJ%csfs)
      end if
      expansion_cfg_JJ%size=0
      expansion_cfg_JJ%nr_of_state=0
!
!
      end subroutine delete_cfg_expansions
!
!-------  function equivalent_JJJK2   ------------
!
      function equivalent_JJJK(csf1, csf2)   result(rez)
!--------------------------------------------------------------------
! This function defines the "equivalency" of
! two csfs in JJ coupling for the formation of
! JK csfs
!(for notion of the "equivalency" see the
! description in the program)
!--------------------------------------------------------------------
      implicit none
      type(csf_LS) :: csf1, csf2
      logical :: rez
      integer::isubc
      rez = .true.
!
      if(csf1%nosubc.lt.2.or.csf2%nosubc.lt.2) then
         write(*,*)'ERROR at subroutine equivalent_JJJK: ',            &
                                  'csf1%nosubc.lt.2.or.csf2%nosubc.lt.2'
         write(*,*)'(you can not use JK coupling with less than 2 ',   &
                                                            'subshelJJ)'
         write(*,*)'program is terminated'
         stop
      end if
      if(csf1%nosubc.eq.csf2%nosubc) then
!
         do isubc=1,csf1%nosubc-1
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            end if
            if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
               rez=.false.
               exit
            end if
            if(.not.csf1%iM1(isubc)==csf2%iM1(isubc)) then
               rez=.false.
               exit
            end if
            if(.not.csf1%iM2(isubc)==csf2%iM2(isubc)) then
               rez=.false.
               exit
            end if
         end do
         isubc=csf1%nosubc
         if(rez) then
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc))        &
                                                            rez=.false.
            if(.not.csf1%subc(isubc)==csf2%subc(isubc)) rez=.false.
            if(.not.csf1%iM2(isubc)==csf2%iM2(isubc)) rez=.false.
         end if
      else
         rez=.false.
      end if
      end function equivalent_JJJK
!
!-------  subroutine form_list_nomach_csfs   ------------
!
      subroutine form_list_nomach_csfs
!--------------------------------------------------------------------
! This subroutine form the list of serial numbers
! in "expansion_cfg_JJ" of "nonequivalent"
! JJ coupling csfs
!(for notion of the "equivalency" see the
! description in the program)
!--------------------------------------------------------------------
      implicit none
      integer:: inomach_counter
      integer, dimension(:), pointer :: temp_list
      integer::icsf_JJ, inew
      logical:: new_csf
!
      inomach_counter=0
      allocate(temp_list(expansion_cfg_JJ%size))
!
      do icsf_JJ=1,expansion_cfg_JJ%size,1
         new_csf = .true.
         do inew = 1, inomach_counter, 1
!write(*,*)'cafs to check: ',icsf_JJ, inew
            if(equivalent_JJJK(expansion_cfg_JJ%csfs(icsf_JJ),         &
              expansion_cfg_JJ%csfs(temp_list(inew)))) new_csf = .false.
!if(equivalent_in_matrix_JJJK(expansion_cfg_JJ%csfs(icsf_JJ),expansion_cfg_JJ%csfs(temp_list(inew)))) new_csf = .false.
         end do
         if(new_csf) then
            inomach_counter = inomach_counter + 1
            if(inomach_counter.gt.expansion_cfg_JJ%size) then
               write(*,*) 'stop at subroutine define_nomach: ',        &
                  'inomach_counter.gt.expansion_cfg_JJ%size'
               stop
            end if
            temp_list(inomach_counter)= icsf_JJ
         end if
      end do !icsf_JJ
!
!write(*,*)' inomach_counter', inomach_counter
!
      nomach_csfs_JJ%list_size = inomach_counter
!
      allocate(nomach_csfs_JJ%items(nomach_csfs_JJ%list_size))
!
      do icsf_JJ=1, nomach_csfs_JJ%list_size, 1
         nomach_csfs_JJ%items(icsf_JJ) = temp_list(icsf_JJ)
      end do !icsf_JJ
      deallocate(temp_list)
!----- print down data -----------
!write(iwrite_log, *)'Nonmaching csfs_JJ:'
!write(iwrite_log, *)'       Nr.   icsf_nr '
!do icsf_JJ=1, nomach_csfs_JJ%list_size, 1
!	 write(iwrite_log, '(4x,i3,2x,i3)') icsf_JJ, nomach_csfs_JJ%items(icsf_JJ)
!end do !icsf_JJ
!write(iwrite_log, *)'  '
      end subroutine form_list_nomach_csfs
!
!***********************************************************************
!                                                                      *
      subroutine form_csfs_JK(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in JK coupling "expansion_cfg_JK",                     *
!     correponding to the one in JJ coupling "expansion_cfg_JJ"        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
      integer:: itype !defines the type of execution
                !itype=1 - to find the number
                !of csfs_JK only,
                !itype.ne.1 - perform full
                !calculation
      integer :: inr_of_csfs_JK   !, icsf_nr
      integer :: icsf_JJ, icsf_JK
      integer :: K, K_min, K_max, J, J_min,J_max
      integer :: isubc_aviable
      integer :: inosubc, isubc
      integer :: iJ_total
      integer :: inr_of_csf_in_expansion_JJ
      isubc_aviable = 2
      if(expansion_cfg_JJ%csfs(1)%nosubc.gt.isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
            expansion_cfg_JJ%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_JK ',    &
            'module transform_JJjj: cfg_Ji_structure%nr_of_subc .gt. ',&
            'isubc_aviable'
         stop
      end if
      iJ_total=states%states(expansion_cfg_JJ%nr_of_state)%J
      call define_number_of_csfs_JK(inr_of_csfs_JK)
      expansion_cfg_JK%size=inr_of_csfs_JK
      if(itype.ne.1) then
         allocate(expansion_cfg_JK%csfs(expansion_cfg_JK%size))
         allocate(expansion_cfg_JK%coeffs(expansion_cfg_JK%size))
         icsf_JK=0
         do icsf_JJ=1, nomach_csfs_JJ%list_size,1
            inr_of_csf_in_expansion_JJ = nomach_csfs_JJ%items(icsf_JJ)
            inosubc = expansion_cfg_JJ%csfs(                           &
               inr_of_csf_in_expansion_JJ)%nosubc
            K_min = abs(expansion_cfg_JJ%csfs(                         &
               inr_of_csf_in_expansion_JJ)%iM2(inosubc-1) -            &
               expansion_cfg_JJ%csfs(inr_of_csf_in_expansion_JJ)%      &
               subc(inosubc)%iL)
            K_max = abs(expansion_cfg_JJ%csfs(                         &
               inr_of_csf_in_expansion_JJ)%iM2(inosubc-1) +            &
               expansion_cfg_JJ%csfs(inr_of_csf_in_expansion_JJ)%      &
               subc(inosubc)%iL)
            do K= K_min, K_max, 2
               J_min = abs(K - expansion_cfg_JJ%csfs(                  &
                  inr_of_csf_in_expansion_JJ)%subc(inosubc)%iS)
               J_max = abs(K + expansion_cfg_JJ%csfs(                  &
                  inr_of_csf_in_expansion_JJ)%subc(inosubc)%iS)
               do J = J_min, J_max, 2
                  if(J.eq.iJ_total) then
                     icsf_JK = icsf_JK + 1
                     if(icsf_JK .gt. expansion_cfg_JK%size) then
                        write(*,*) 'stop at subroutine form_csfs_JK: ',&
                           'icsf_nr.gt.expansion_cfg_JK%size'
                        stop
                     end if
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc=inosubc
                     allocate(expansion_cfg_JK%csfs(icsf_JK)%subc_cfg( &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc))
                     allocate(expansion_cfg_JK%csfs(icsf_JK)%subc(     &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc))
                     allocate(expansion_cfg_JK%csfs(icsf_JK)%iM1(      &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc))
                     allocate(expansion_cfg_JK%csfs(icsf_JK)%iM2(      &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc))
                     allocate(expansion_cfg_JK%csfs(icsf_JK)%iJ(       &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc))
                     expansion_cfg_JK%csfs(icsf_JK)%iJ(                &
                     expansion_cfg_JK%csfs(icsf_JK)%nosubc) = 0
                     do isubc=1,expansion_cfg_JK%csfs(icsf_JK)%        &
                                                             nosubc-1,1
                         expansion_cfg_JK%csfs(icsf_JK)%iJ(isubc) = 0
                        expansion_cfg_JK%csfs(icsf_JK)%subc_cfg(isubc)=&
                           expansion_cfg_JJ%csfs(                      &
                           inr_of_csf_in_expansion_JJ)%subc_cfg(isubc)
                        expansion_cfg_JK%csfs(icsf_JK)%subc(isubc) =   &
                        expansion_cfg_JJ%csfs(                         &
                        inr_of_csf_in_expansion_JJ)%subc(isubc)
                        expansion_cfg_JK%csfs(icsf_JK)%iM1(isubc) =    &
                        expansion_cfg_JJ%csfs(                         &
                        inr_of_csf_in_expansion_JJ)%iM1(isubc)
                        expansion_cfg_JK%csfs(icsf_JK)%iM2(isubc) =    &
                        expansion_cfg_JJ%csfs(                         &
                        inr_of_csf_in_expansion_JJ)%iM2(isubc)
                     end do
                     isubc = expansion_cfg_JK%csfs(icsf_JK)%nosubc
                     expansion_cfg_JK%csfs(icsf_JK)%subc_cfg(isubc) =  &
                     expansion_cfg_JJ%csfs(inr_of_csf_in_expansion_JJ)%&
                     subc_cfg(isubc)
                     expansion_cfg_JK%csfs(icsf_JK)%subc(isubc) =      &
                     expansion_cfg_JJ%csfs(inr_of_csf_in_expansion_JJ)%&
                     subc(isubc)
                     expansion_cfg_JK%csfs(icsf_JK)%iM1(isubc)= K
                     expansion_cfg_JK%csfs(icsf_JK)%iM2(isubc) =       &
                     expansion_cfg_JJ%csfs(inr_of_csf_in_expansion_JJ)%&
                     iM2(isubc)
                  end if
               end do
            end do
         end do
      end if
      contains
!
!--- 	 subroutine  define_number_of_csfs_JK  ---------
!
         subroutine  define_number_of_csfs_JK(irez)
!--------------------------------------------------------------------
! This subroutine defines the number of csfs in
! JK coupling
!--------------------------------------------------------------------
         implicit none
         integer, intent(out)::irez
         integer::K, K_min, K_max, J, J_min,J_max
         integer::isubc_aviable
         integer::icsf
         integer :: icsf_nr_in_expansion_cfg_JJ, inosubc
         isubc_aviable = 2
         if(expansion_cfg_JJ%csfs(1)%nosubc .gt. isubc_aviable) then
            write(*,*) 'STOP at subroutine define_number_of_csfs_JK ', &
              'module transform_JJjj: cfg_Ji_structure%nr_of_subc.gt.',&
              'isubc_aviable'
            stop
         end if
         irez=0
         do icsf=1, nomach_csfs_JJ%list_size,1
            icsf_nr_in_expansion_cfg_JJ = nomach_csfs_JJ%items(icsf)
            inosubc = expansion_cfg_JJ%csfs(                           &
               icsf_nr_in_expansion_cfg_JJ)%nosubc
            K_min = abs(expansion_cfg_JJ%csfs(                         &
               icsf_nr_in_expansion_cfg_JJ)%iM2(inosubc-1) -           &
               expansion_cfg_JJ%csfs(icsf_nr_in_expansion_cfg_JJ)%     &
               subc(inosubc)%iL)
            K_max = abs(expansion_cfg_JJ%csfs(                         &
               icsf_nr_in_expansion_cfg_JJ)%iM2(inosubc-1) +           &
               expansion_cfg_JJ%csfs(icsf_nr_in_expansion_cfg_JJ)%     &
               subc(inosubc)%iL)
            do K= K_min, K_max, 2
               J_min=abs(K - expansion_cfg_JJ%csfs(                    &
                  icsf_nr_in_expansion_cfg_JJ)%subc(inosubc)%iS)
               J_max=abs(K + expansion_cfg_JJ%csfs(                    &
                  icsf_nr_in_expansion_cfg_JJ)%subc(inosubc)%iS)
               do J = J_min, J_max, 2
                  if(J.eq.iJ_total) irez = irez + 1
               end do
            end do
         end do
         end subroutine  define_number_of_csfs_JK
      end subroutine form_csfs_JK
!
!--- 	 subroutine  matrix_JJ_JK  ---------
!
      subroutine matrix_JJ_JK(icsf_JJ,icsf_JK,rez)
!--------------------------------------------------------------------
! This subroutine calculates the transformation
! matrix element between two csfs
! (in JJ and JK couplings)
!--------------------------------------------------------------------
      implicit none
      integer,intent(in):: icsf_JJ,icsf_JK
      real(kind=dp),intent(out)::rez
      integer::K,iSi,J,iJi,iJ_im1,iLi
      integer::nosubc !number of subshells
      rez = ZERO_dp
      if(equivalent_JJJK(expansion_cfg_JJ%csfs(icsf_JJ),               &
                                  expansion_cfg_JK%csfs(icsf_JK))) then
         nosubc=expansion_cfg_JJ%csfs(icsf_JJ)%nosubc
         K=expansion_cfg_JK%csfs(icsf_JK)%iM1(nosubc)
         iSi=expansion_cfg_JJ%csfs(icsf_JJ)%subc(nosubc)%iS
         J=states%states(expansion_cfg_JJ%nr_of_state)%J
         iJi=expansion_cfg_JJ%csfs(icsf_JJ)%iM1(nosubc)
         iJ_im1=expansion_cfg_JJ%csfs(icsf_JJ)%iM2(nosubc-1)
         iLi=expansion_cfg_JJ%csfs(icsf_JJ)%subc(nosubc)%iL
         call JJJKPER(K,iSi,J,iJi,iJ_im1,iLi,rez)
      end if
      end subroutine matrix_JJ_JK
!
!--- 	 subroutine  dotransform_JJjj  ---------
!
      subroutine dotransform_JJJK
!--------------------------------------------------------------------
! This subroutine calculates the weights of the
! expansions in the JK coupling
!--------------------------------------------------------------------
      implicit none
      real(kind=dp)::coeff_JK,coeff_JJ, melement, sum_of_state
      integer::icsf_JJ,icsf_JK
      sum_of_state = 0.
      do icsf_JK=1,expansion_cfg_JK%size
         coeff_JK = 0.
         do icsf_JJ=1, expansion_cfg_JJ%size
            call matrix_JJ_JK(icsf_JJ,icsf_JK,melement)
            coeff_JJ=expansion_cfg_JJ%coeffs(icsf_JJ)
            coeff_JK = coeff_JK + coeff_JJ*melement
         end do
         if((dabs(coeff_JK)-dp_coeff_precision).gt.ONE_dp) then
            write(*,*)'possible error at subroutine dotransform_JJJK:',&
                                                  ' coeff_JK=',coeff_JK
            write(iwrite_log,*)'possible error at subroutine ',        &
                                 'dotransform_JJJK: coeff_JK=',coeff_JK
         end if
         expansion_cfg_JK%coeffs(icsf_JK)= coeff_JK
         sum_of_state = sum_of_state + coeff_JK*coeff_JK
      end do
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_JJJK: ',  &
                                           'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
                         'dotransform_JJJK: sum_of_state=',sum_of_state
      end if
      end subroutine dotransform_JJJK
!
!------   subroutine print_cfg_JJJK      -------------------
!
      subroutine print_cfg_JJJK(itype)
!--------------------------------------------------------------------
! This subroutine prints the oneconfigurational
! expansions to the unit "iwrite_cfg_expansions_JK"
! (see module "Coupling_constants")
!--------------------------------------------------------------------
      implicit none
      integer, intent(in) ::itype
      integer             :: icsf, isubc, j
      character(len=1) :: CVAL
      character(len=4) :: JVAL
!
!     JJ - Coupling
!
      if(itype.gt.0) then
         write(iwrite_cfg_expansions_JK,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JK,*) &
         'state Nr.',expansion_cfg_JJ%nr_of_state
         write(iwrite_cfg_expansions_JK,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_JJ%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_JJ%nr_of_state)%energy
         write(iwrite_cfg_expansions_JK,'(3x,a29,i2)')               &
         'expansion size (JJ coupling): ', expansion_cfg_JJ%size
         write(iwrite_cfg_expansions_JK,*)''
         if(associated(expansion_cfg_JJ%csfs)) then
            do icsf=1,expansion_cfg_JJ%size
               if(associated(expansion_cfg_JJ%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_JK,*)'csf Nr.'
                  if(expansion_cfg_JJ%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK,                  &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf, &
                     expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%in,      &
                   CVAL(1,expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%iN_big,  &
                     'coeff:',expansion_cfg_JJ%coeffs(icsf)
                  else if (expansion_cfg_JJ%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK,                   &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%il),&
                    expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%iN_big,    &
                    j=1,2),                                            &
                    'coeff:',expansion_cfg_JJ%coeffs(icsf)
                  else
                    STOP &
                  'To many couled shells in Coupling_transform_cfg_JJJK'
                  end if
               else
                  write(iwrite_cfg_expansions_JK,'(9x,a51)') &
                  'expansion_cfg_JJ%csfs(icsf)%subc_cfg NOT associated'
               end if
                  if(associated(expansion_cfg_JJ%csfs(icsf)%subc)   &
                  .and.associated(expansion_cfg_JJ%csfs(icsf)%iM1)  &
                  .and.associated(expansion_cfg_JJ%csfs(icsf)%iM2)) then
                  if(expansion_cfg_JJ%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_JJ%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_JJ%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_JJ%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_JJ%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK,                &
                     '(12x,2(1x,i2,a1,i1),3x,"[",a4,","a4,"]",a4)') &
                     (expansion_cfg_JJ%csfs(icsf)%subc(j)%iS+1,     &
                     CVAL(2,expansion_cfg_JJ%csfs(icsf)%subc(j)%iL),&
                     expansion_cfg_JJ%csfs(icsf)%subc(j)%inr,       &
                     j=1,2),                                        &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM2(1)),      &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM1(2)),      &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM2(2))
                  end if
               else
                  write(iwrite_cfg_expansions_JK,'(9x,a63)')  &
       'expansion_cfg_JJ%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JK,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JK,*)  &
            'expansion_cfg_JJ%csfs NOT associated'
         end if
      end if
!
!     JK - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_JK,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JK,*) &
         'state Nr.',expansion_cfg_JK%nr_of_state
         write(iwrite_cfg_expansions_JK,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_JK%nr_of_state)%J), &
         ' Energy =', states%states(expansion_cfg_JK%nr_of_state)%energy
         write(iwrite_cfg_expansions_JK,'(3x,a29,i2)')               &
         'expansion size (JK coupling): ', expansion_cfg_JK%size
         write(iwrite_cfg_expansions_JK,*)''
         if(associated(expansion_cfg_JK%csfs)) then
            do icsf=1,expansion_cfg_JK%size
               if(associated(expansion_cfg_JK%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_JK,*)'csf Nr.'
                  if(expansion_cfg_JK%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK,                  &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf, &
                     expansion_cfg_JK%csfs(icsf)%subc_cfg(1)%in,      &
                   CVAL(1,expansion_cfg_JK%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_JK%csfs(icsf)%subc_cfg(1)%iN_big,  &
                     'coeff:',expansion_cfg_JK%coeffs(icsf)
                  else if (expansion_cfg_JK%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK,                   &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_JK%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_JK%csfs(icsf)%subc_cfg(j)%il),&
                    expansion_cfg_JK%csfs(icsf)%subc_cfg(j)%iN_big,    &
                    j=1,2),                                            &
                    'coeff:',expansion_cfg_JK%coeffs(icsf)
                  else
                    STOP &
                  'To many couled shells in Coupling_transform_cfg_JJJK'
                  end if
               else
                  write(iwrite_cfg_expansions_JK,'(9x,a51)') &
                  'expansion_cfg_JK%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_JK%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_JK%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_JK%csfs(icsf)%iM2)) then
                  if(expansion_cfg_JK%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_JK%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_JK%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_JK%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_JK%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK,                &
                     '(12x,2(1x,i2,a1,i1),3x,a4," ["a4,"]",a4)') &
                     (expansion_cfg_JK%csfs(icsf)%subc(j)%iS+1,     &
                     CVAL(2,expansion_cfg_JK%csfs(icsf)%subc(j)%iL),&
                     expansion_cfg_JK%csfs(icsf)%subc(j)%inr,       &
                     j=1,2),                                        &
                     JVAL(expansion_cfg_JK%csfs(icsf)%iM2(1)),      &
                     JVAL(expansion_cfg_JK%csfs(icsf)%iM1(2)),      &
                     JVAL(expansion_cfg_JK%csfs(icsf)%iM2(2))
                  end if
               else
                  write(iwrite_cfg_expansions_JK,'(9x,a63)') &
       'expansion_cfg_JK%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JK,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JK,*) &
            'expansion_cfg_JK%csfs NOT associated'
         end if
         write(iwrite_cfg_expansions_JK,*) &
         '-------------------------------------'
      end if
      end subroutine print_cfg_JJJK
!
      end module Coupling_transform_cfg_JJJK
