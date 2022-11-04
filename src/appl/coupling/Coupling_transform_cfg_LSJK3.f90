!
!***********************************************************************
!
      module Coupling_transform_cfg_LSJK3
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
      public  :: main_cfg_lsjk3
      public  :: count_nr_of_csfs_JK3
      public  :: delete_cfg_expansions
      private :: form_list_nomach_csfs
      private :: form_csfs_JK3
      private :: matrix_LS_JK3
      private :: dotransform_LSJK3
      private :: print_cfg_LSJK3
      private :: equivalent_LSJK3
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
      type(expansion), public :: expansion_cfg_JK3
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_cfg_lsjk3(print_level)
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
      expansion_cfg_JK3%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=0
      call form_csfs_JK3(itype)
      call dotransform_LSJK3
      if(print_level.gt.1) call print_cfg_LSJK3(2)
      deallocate(nomach_csfs_ls%items, STAT = error)
      end subroutine main_cfg_lsjk3
!
!***********************************************************************
!                                                                      *
      subroutine count_nr_of_csfs_JK3(nr_of_csfs)
!                                                                      *
!     This subroutine counts the number of oneconfigurational          *
!     expansion in JK3 coupling                                        *
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
      if(expansion_cfg_LS%csfs(1)%nosubc /= 3) then
         write(*,*)                                                    &
            'ERROR at subroutine count_nr_of_csfs_JK3: nosubc /= 3'
         write(*,*)                                                    &
            '(you can use JK3 coupling with subshells number  = 3)'
         write(*,*) 'At presint number of subshells is',               &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*)'program is terminated'
         stop
      end if
      expansion_cfg_JK3%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=1
      call form_csfs_JK3(itype)
      nr_of_csfs = expansion_cfg_JK3%size
      if(associated(nomach_csfs_ls%items))                             &
                             deallocate(nomach_csfs_ls%items,STAT=error)
      end subroutine count_nr_of_csfs_JK3
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_expansions
!                                                                      *
!     This subroutine deallocates the arrays of                        *
!     oneconfigurational expansions "expansion_cfg_JK3" and           *
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
!delete expansion_cfg_JK3
      if(associated(expansion_cfg_JK3%coeffs))                        &
                                     deallocate(expansion_cfg_JK3%csfs)
      if(associated(expansion_cfg_JK3%csfs)) then
         do icsf=1, expansion_cfg_JK3%size
            if(associated(expansion_cfg_JK3%csfs(icsf)%subc_cfg))     &
                      deallocate(expansion_cfg_JK3%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_JK3%csfs(icsf)%subc))         &
                          deallocate(expansion_cfg_JK3%csfs(icsf)%subc)
            if(associated(expansion_cfg_JK3%csfs(icsf)%iM1))          &
                           deallocate(expansion_cfg_JK3%csfs(icsf)%iM1)
            if(associated(expansion_cfg_JK3%csfs(icsf)%iM2))          &
                           deallocate(expansion_cfg_JK3%csfs(icsf)%iM2)
            if(associated(expansion_cfg_JK3%csfs(icsf)%iJ))           &
                           deallocate(expansion_cfg_JK3%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_JK3%csfs)
      end if
      expansion_cfg_JK3%size=0
      expansion_cfg_JK3%nr_of_state=0
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
      end subroutine delete_cfg_expansions
!
!***********************************************************************
!                                                                      *
      function equivalent_LSJK3(csf1, csf2)                result(rez)
!                                                                      *
!     This function defines the "equivalency" of                       *
!     two csfs in LS coupling for the formation of JK3 csfs            *
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
      if(csf1%nosubc /= 3 .or. csf2%nosubc /= 3) then
         write(*,*)'ERROR at subroutine equivalent_LSJK3: ',          &
                   'csf1%nosubc /= 3 .or. csf2%nosubc /= 3'
         write(*,*)'(you can use JK3 coupling with shell number 3)'
         write(*,*)'program will be terminated'
         stop
      end if
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
               rez=.false.
               exit
            end if
!            if(isubc /= 2) then
            if(isubc == 1) then
               if(.not.csf1%iM1(isubc)==csf2%iM1(isubc)) then
                  rez=.false.
                  exit
               else if(.not.csf1%iM2(isubc)==csf2%iM2(isubc)) then
                  rez=.false.
                  exit
               end if
            end if
         end do
      else
         rez=.false.
      end if
      end function equivalent_LSJK3
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
            if(equivalent_LSJK3(expansion_cfg_LS%csfs(icsf_LS),       &
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
      subroutine form_csfs_JK3(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in JK3 coupling "expansion_cfg_JK3",                 *
!     correponding to the one in LS coupling "expansion_cfg_LS"        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
!type(cfg_Ji_lists),intent(in)::cfg_Ji_structure
      integer:: itype !defines the type of execution
                      !itype=1 - to find the number
                      !of csfs_JK3 only,
                      !itype.ne.1 - perform full
                      !calculation
      integer :: inr_of_csfs_JK3
      integer :: icsf_LS, icsf_JK3
      integer :: L, L_min, L_max, IS, IS_min, IS_max
      integer :: J_12, J_12_min, J_12_max, K, K_min, K_max
      integer :: isubc_aviable, ITTK
      integer :: inosubc, isubc
      integer :: iJ_total
      integer :: icsf_nr
      isubc_aviable = 3
      if(expansion_cfg_LS%csfs(1)%nosubc.gt.isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_JK3 ',  &
            'module transform_lsjk3: cfg_Ji_structure%nr_of_subc.gt.',&
            'isubc_aviable'
         stop
      end if
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
      call define_number_of_csfs_JK3(inr_of_csfs_JK3)
      expansion_cfg_JK3%size=inr_of_csfs_JK3
      if(itype.ne.1) then
         allocate(expansion_cfg_JK3%csfs(expansion_cfg_JK3%size))
         allocate(expansion_cfg_JK3%coeffs(expansion_cfg_JK3%size))
         icsf_JK3=0
         do icsf_LS=1, nomach_csfs_ls%list_size,1
           icsf_nr = nomach_csfs_ls%items(icsf_LS)
           inosubc = expansion_cfg_LS%csfs(icsf_nr)%nosubc
           L_min = iabs(                                               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL-               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iL)
           L_max =                                                     &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL+               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iL
           IS_min = iabs(                                              &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS-               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iS)
           IS_max =                                                    &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS+               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iS
           J_12_min = iabs(                                            &
              expansion_cfg_LS%csfs(icsf_nr)%iM1(1)-                   &
              expansion_cfg_LS%csfs(icsf_nr)%iM2(1))
           J_12_max =                                                  &
              expansion_cfg_LS%csfs(icsf_nr)%iM1(1)+                   &
              expansion_cfg_LS%csfs(icsf_nr)%iM2(1)
           do L = L_min, L_max, 2
              do J_12 = J_12_min, J_12_max, 2
                 K_min = iabs(L - J_12)
                 K_max = L + J_12
                 do K = K_min, K_max, 2
                    do IS = IS_min, IS_max, 2
                       if(ITTK(2*K,2*IS,2*iJ_total) == 1) then
                          icsf_JK3 = icsf_JK3 + 1
                          if(icsf_JK3.gt.expansion_cfg_JK3%size) then
                             write(8,*)                                &
                               'stop at subroutine form_csfs_JK3:',   &
                               ' icsf_nr.gt.expansion_cfg_JK3%size'
                             stop
                          end if
                          expansion_cfg_JK3%csfs(icsf_JK3)%nosubc =  &
                                                                inosubc
                          allocate(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          subc_cfg(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          nosubc))
                          allocate(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          subc(expansion_cfg_JK3%csfs(icsf_JK3)%     &
                          nosubc))
                          allocate(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          iM1(expansion_cfg_JK3%csfs(icsf_JK3)%      &
                                                               nosubc))
                          allocate(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          iM2(expansion_cfg_JK3%csfs(icsf_JK3)%      &
                                                                nosubc))
                          allocate(expansion_cfg_JK3%csfs(icsf_JK3)% &
                          iJ(expansion_cfg_JK3%csfs(icsf_JK3)%nosubc))
                          do isubc=1,expansion_cfg_JK3%csfs(          &
                                                     icsf_JK3)%nosubc,1
                             expansion_cfg_JK3%csfs(icsf_JK3)%       &
                             subc_cfg(isubc) = expansion_cfg_LS%csfs(  &
                             icsf_nr)%subc_cfg(isubc)
                             expansion_cfg_JK3%csfs(icsf_JK3)%       &
                             subc(isubc) = expansion_cfg_LS%csfs(      &
                             icsf_nr)%subc(isubc)
                             if(isubc == 1) then
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                               iM1(isubc) = expansion_cfg_LS%csfs(     &
                               icsf_nr)%iM1(isubc)
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                               iM2(isubc) = expansion_cfg_LS%csfs(     &
                               icsf_nr)%iM2(isubc)
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                                                        iJ(isubc)=J_12
                             else if(isubc == 2) then
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                                                         iM1(isubc)=L
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                                                        iM2(isubc)=IS
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                                                        iJ(isubc)= K
                             else if(isubc == 3) then
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                               iM1(isubc) = expansion_cfg_LS%csfs(     &
                               icsf_nr)%iM1(isubc)
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                               iM2(isubc) = expansion_cfg_LS%csfs(     &
                               icsf_nr)%iM2(isubc)
                               expansion_cfg_JK3%csfs(icsf_JK3)%     &
                                                      iJ(isubc)=iJ_total
                             end if
                          end do
                       end if
                    end do
                 end do
              end do
           end do
         end do
!
      end if
!
!write(*,*) '      end subroutine form_csfs_JK3'
!write(iwrite_log,*) '      end subroutine form_csfs_JK3'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine  define_number_of_csfs_JK3(irez)
!                                                                      *
!     This subroutine defines the number of csfs in JK3 coupling       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
         implicit none
         integer, intent(out) :: irez
         integer :: L, L_min, L_max, IS, IS_min, IS_max
         integer :: J_12, J_12_min, J_12_max, K, K_min, K_max
         integer :: isubc_aviable
         integer :: icsf, icsf_nr, ITTK
         isubc_aviable = 3
         if(expansion_cfg_LS%csfs(1)%nosubc .gt. isubc_aviable) then
           write(*,*) 'STOP at subroutine define_number_of_csfs_JK3 ',&
             'module transform_lslsJ3:cfg_Ji_structure%nr_of_subc.gt.',&
              'isubc_aviable'
           stop
         end if
         irez=0
         do icsf=1, nomach_csfs_ls%list_size,1
           icsf_nr = nomach_csfs_ls%items(icsf)
           L_min = iabs(                                               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL-               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iL)
           L_max =                                                     &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL+               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iL
           IS_min = iabs(                                              &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS-               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iS)
           IS_max =                                                    &
              expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS+               &
              expansion_cfg_LS%csfs(icsf_nr)%subc(3)%iS
           J_12_min = iabs(                                            &
              expansion_cfg_LS%csfs(icsf_nr)%iM1(1)-                   &
              expansion_cfg_LS%csfs(icsf_nr)%iM2(1))
           J_12_max =                                                  &
              expansion_cfg_LS%csfs(icsf_nr)%iM1(1)+                   &
              expansion_cfg_LS%csfs(icsf_nr)%iM2(1)
           do L = L_min, L_max, 2
              do J_12 = J_12_min, J_12_max, 2
                 K_min = iabs(L - J_12)
                 K_max = L + J_12
                 do K = K_min, K_max, 2
                    do IS = IS_min, IS_max, 2
                       if(ITTK(2*K,2*IS,2*iJ_total) == 1)         &
                                                       irez = irez + 1
                    end do
                 end do
              end do
           end do
         end do
         end subroutine  define_number_of_csfs_JK3
      end subroutine form_csfs_JK3 !
!***********************************************************************
!                                                                      *
      subroutine matrix_LS_JK3(icsf_LS,icsf_JK3,rez)
!                                                                      *
!     This subroutine calculates the transformation                    *
!     matrix element between two csfs  (in LS and JK3 couplings)       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer,intent(in) :: icsf_LS, icsf_JK3
      real(kind=dp),intent(out) :: rez
      integer :: L_12, L_3, L_123, L_4, L_34, L
      integer :: IS_12, IS_3, IS_123, IS_4, IS_34, IS
      integer :: J_12, K, J
      integer :: nosubc
      nosubc=expansion_cfg_LS%csfs(1)%nosubc
      rez = ZERO_dp
      if(equivalent_LSJK3(expansion_cfg_LS%csfs(icsf_LS),             &
         expansion_cfg_JK3%csfs(icsf_JK3))) then
         L_12   = expansion_cfg_LS%csfs(icsf_LS)%iM1(1)
         L_3    = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iL
         L_123  = expansion_cfg_LS%csfs(icsf_LS)%iM1(2)
         L_4    = expansion_cfg_LS%csfs(icsf_LS)%subc(3)%iL
         L_34   = expansion_cfg_JK3%csfs(icsf_JK3)%iM1(2)
         L      = expansion_cfg_LS%csfs(icsf_LS)%iM1(3)
         IS_12  = expansion_cfg_LS%csfs(icsf_LS)%iM2(1)
         IS_3   = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iS
         IS_123 = expansion_cfg_LS%csfs(icsf_LS)%iM2(2)
         IS_4   = expansion_cfg_LS%csfs(icsf_LS)%subc(3)%iS
         IS_34  = expansion_cfg_JK3%csfs(icsf_JK3)%iM2(2)
         IS     = expansion_cfg_LS%csfs(icsf_LS)%iM2(3)
         J_12   = expansion_cfg_JK3%csfs(icsf_JK3)%iJ(1)
         K      = expansion_cfg_JK3%csfs(icsf_JK3)%iJ(2)
         J      = expansion_cfg_JK3%csfs(icsf_JK3)%iJ(3)
         call JK3PER(L_12, L_3, L_123, L_4, L_34, L,                  &
                     IS_12,IS_3,IS_123,IS_4,IS_34,IS,J_12,K,J,rez)
      end if
      end subroutine matrix_LS_JK3
!
!***********************************************************************
!                                                                      *
      subroutine dotransform_LSJK3
!                                                                      *
!     This subroutine calculates the weights of the                    *
!     expansions in the JK3 coupling                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp)::coeff_JK3,coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_JK3
!
      melement=0.
      sum_of_state = 0.
      do icsf_JK3=1,expansion_cfg_JK3%size
         coeff_JK3 = 0.
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_JK3(icsf_LS,icsf_JK3,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_JK3 = coeff_JK3 + coeff_LS*melement
         end do
         if((dabs(coeff_JK3)-dp_coeff_precision).gt.ONE_dp) then
          write(*,*)'possible error at subroutine dotransform_LSJK3:',&
               ' coeff_JK3=',coeff_JK3
            write(iwrite_log,*)'possible error at subroutine ',        &
               'dotransform_LSJK3: coeff_JK3=',coeff_JK3
         end if
         expansion_cfg_JK3%coeffs(icsf_JK3)= coeff_JK3
         sum_of_state = sum_of_state + coeff_JK3*coeff_JK3
      end do
!
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_LSJK3: ',&
            'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
            'dotransform_LSJK3: sum_of_state=',sum_of_state
      end if
!
      end subroutine dotransform_LSJK3
!
!***********************************************************************
!                                                                      *
      subroutine print_cfg_LSJK3(itype)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the unit "iwrite_cfg_expansions_JK3"              *
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
         write(iwrite_cfg_expansions_JK3,*)  &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JK3,*)  &
         'state Nr.',expansion_cfg_LS%nr_of_state
         write(iwrite_cfg_expansions_JK3,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_LS%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LS%nr_of_state)%energy
         write(iwrite_cfg_expansions_JK3,'(3x,a30,i2)')               &
         'expansion size (LS coupling): ', expansion_cfg_LS%size
         write(iwrite_cfg_expansions_JK3,*)''
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1,expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_JK3,*)      &
                     'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf,  &
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%in,       &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%iN_big,   &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(1x,i5,3x,3(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,3),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
                 'To many couled shells in Coupling_transform_cfg_LSJK3'
                  end if
               else
                  write(iwrite_cfg_expansions_JK3,'(9x,a51)')   &
                  'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LS%csfs(icsf)%subc)  &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM1) &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JK3,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_LS%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LS%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(12x,7(1x,i2,a1,i1))')                         &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,      &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL),&
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,       &
                      j=1,2),                                        &
                      expansion_cfg_LS%csfs(icsf)%iM2(2)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(2))
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(12x,3(1x,i2,a1,i1),3(1x,i2,a1))')              &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,       &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL), &
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,        &
                      j=1,3),                                         &
                      (expansion_cfg_LS%csfs(icsf)%iM2(j)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(j)),j=2,3)
                  end if
               else
                   write(iwrite_cfg_expansions_JK3,'(9x,a63)') &
       'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JK3,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JK3,*) &
            'expansion_cfg_LS%csfs NOT associated'
         end if
      end if
!
!     JK3 - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_JK3,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JK3,*) &
         'state Nr.',expansion_cfg_JK3%nr_of_state
         write(iwrite_cfg_expansions_JK3,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_JK3%nr_of_state)%J), &
        ' Energy =',states%states(expansion_cfg_JK3%nr_of_state)%energy
         write(iwrite_cfg_expansions_JK3,'(3x,a32,i2)')               &
         'expansion size (JK3 coupling): ', expansion_cfg_JK3%size
         write(iwrite_cfg_expansions_JK3,*)''
         if(associated(expansion_cfg_JK3%csfs)) then
            do icsf=1,expansion_cfg_JK3%size
               if(associated(expansion_cfg_JK3%csfs(icsf)%subc_cfg)) then
!                  if(dabs(expansion_cfg_JK3%coeffs(icsf)) <=          &
!                                              dp_coeff_precision) cycle
                  if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(1x,i5,3x,3(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_JK3%csfs(icsf)%subc_cfg(j)%in,    &
                  CVAL(1,expansion_cfg_JK3%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_JK3%csfs(icsf)%subc_cfg(j)%iN_big, &
                     j=1,3),                                           &
                     'coeff:',expansion_cfg_JK3%coeffs(icsf)
                  end if
               else
                  write(iwrite_cfg_expansions_JK3,'(9x,a51)')         &
                 'expansion_cfg_JK3%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_JK3%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_JK3%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_JK3%csfs(icsf)%iM2)) then
                  if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_JK3,                 &
                     '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,i2,a1,i1,")", &
                     i2,a1,1x,a4,";",1x,a4," [",a4," ]")')           &
                     (expansion_cfg_JK3%csfs(icsf)%subc(j)%iS+1,      &
                     CVAL(2,expansion_cfg_JK3%csfs(icsf)%subc(j)%iL), &
                     expansion_cfg_JK3%csfs(icsf)%subc(j)%inr,        &
                     j=1,3),                                          &
                     expansion_cfg_JK3%csfs(icsf)%iM2(2)+1,           &
                     CVAL(2,expansion_cfg_JK3%csfs(icsf)%iM1(2)),     &
                     JVAL(expansion_cfg_JK3%csfs(icsf)%iJ(1)),        &
                     JVAL(expansion_cfg_JK3%csfs(icsf)%iJ(2)),        &
                     JVAL(states%states(expansion_cfg_LS%nr_of_state)%J)
                  end if
               else
                  write(iwrite_cfg_expansions_JK3,'(9x,a63)')         &
                     'expansion_cfg_JK3%csfs(icsf)%subc (',           &
                     'or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JK3,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JK3,*) &
            'expansion_cfg_JK3%csfs NOT associated'
         end if
            write(iwrite_cfg_expansions_JK3,*) &
            '-------------------------------------'
      end if
      end subroutine print_cfg_LSJK3
!
      end module Coupling_transform_cfg_LSJK3
