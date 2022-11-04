!
!***********************************************************************
!
      module Coupling_transform_cfg_LScLSJ3
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
      use Coupling_inside_shell
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: main_cfg_lsclsJ3
      public  :: count_nr_of_csfs_cLSJ3
      public  :: delete_cfg_expansions
      private :: form_list_nomach_csfs
      private :: form_csfs_cLSJ3
      private :: matrix_LS_cLSJ3
      private :: dotransform_LScLSJ3
      private :: print_cfg_LScLSJ3
      private :: equivalent_LScLSJ3
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
      type(expansion), public :: expansion_cfg_cLSJ3
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_cfg_lsclsJ3(print_level)
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
      expansion_cfg_cLSJ3%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=0
      call form_csfs_cLSJ3(itype)
      call dotransform_LScLSJ3
      if(print_level.gt.1) call print_cfg_LScLSJ3(2)
      deallocate(nomach_csfs_ls%items, STAT = error)
      end subroutine main_cfg_lsclsJ3
!
!***********************************************************************
!                                                                      *
      subroutine count_nr_of_csfs_cLSJ3(nr_of_csfs)
!                                                                      *
!     This subroutine counts the number of oneconfigurational          *
!     expansion in cLSJ3 coupling                                      *
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
            'ERROR at subroutine count_nr_of_csfs_cLSJ3: nosubc /= 3'
         write(*,*)                                                    &
            '(you can use cLSJ3 coupling with subshells number  = 3)'
         write(*,*) 'At presint number of subshells is',               &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*)'program is terminated'
         stop
      end if
      expansion_cfg_cLSJ3%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=1
      call form_csfs_cLSJ3(itype)
      nr_of_csfs = expansion_cfg_cLSJ3%size
      if(associated(nomach_csfs_ls%items))                             &
                             deallocate(nomach_csfs_ls%items,STAT=error)
      end subroutine count_nr_of_csfs_cLSJ3
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_expansions
!                                                                      *
!     This subroutine deallocates the arrays of                        *
!     oneconfigurational expansions "expansion_cfg_cLSJ3" and          *
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
!delete expansion_cfg_cLSJ3
      if(associated(expansion_cfg_cLSJ3%coeffs))                       &
                                    deallocate(expansion_cfg_cLSJ3%csfs)
      if(associated(expansion_cfg_cLSJ3%csfs)) then
         do icsf=1, expansion_cfg_cLSJ3%size
            if(associated(expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg))    &
                     deallocate(expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_cLSJ3%csfs(icsf)%subc))        &
                         deallocate(expansion_cfg_cLSJ3%csfs(icsf)%subc)
            if(associated(expansion_cfg_cLSJ3%csfs(icsf)%iM1))         &
                          deallocate(expansion_cfg_cLSJ3%csfs(icsf)%iM1)
            if(associated(expansion_cfg_cLSJ3%csfs(icsf)%iM2))         &
                          deallocate(expansion_cfg_cLSJ3%csfs(icsf)%iM2)
            if(associated(expansion_cfg_cLSJ3%csfs(icsf)%iJ))          &
                          deallocate(expansion_cfg_cLSJ3%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_cLSJ3%csfs)
      end if
      expansion_cfg_cLSJ3%size=0
      expansion_cfg_cLSJ3%nr_of_state=0
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
      function equivalent_LScLSJ3(csf1, csf2)               result(rez)
!                                                                      *
!     This function defines the "equivalency" of                       *
!     two csfs in LS coupling for the formation of cLSJ3 csfs          *
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
         write(*,*)'ERROR at subroutine equivalent_LScLSJ3: ',         &
                   'csf1%nosubc /= 3 .or. csf2%nosubc /= 3'
         write(*,*)'(you can use cLSJ3 coupling with shell number 3)'
         write(*,*)'program will be terminated'
         write(*,*)'csf1%nosubc, csf2%nosubc',csf1%nosubc, csf2%nosubc
         stop
      end if
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
               if(isubc /= 1) rez=.false.
               exit
            end if
         end do
      else
         rez=.false.
      end if
      end function equivalent_LScLSJ3
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
            if(equivalent_LScLSJ3(expansion_cfg_LS%csfs(icsf_LS),      &
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
!      write(99, *)'Nonmaching csfs_LS:'
!      write(99, *)'       Nr.   icsf_nr '
!      do icsf_LS=1, nomach_csfs_LS%list_size, 1
!	 write(99, '(4x,i3,2x,i3)') icsf_LS, nomach_csfs_LS%items(icsf_LS)
!      end do  !icsf_LS
!      write(99, *)'  '
!
      end subroutine form_list_nomach_csfs
!
!***********************************************************************
!                                                                      *
      subroutine form_csfs_cLSJ3(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in cLSJ3 coupling "expansion_cfg_cLSJ3",               *
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
                      !of csfs_cLSJ3 only,
                      !itype.ne.1 - perform full
                      !calculation
!GG!      real(kind=dp) :: coef
      integer :: inr_of_csfs_cLSJ3
      integer :: icsf_LS, icsf_cLSJ3
      integer :: L, L_min, L_max, IS, IS_min, IS_max
      integer :: J_1, J_1_min, J_1_max
      integer :: J_12, J_12_min, J_12_max
      integer :: isubc_aviable, ITTK
      integer :: inosubc, isubc
      integer :: iJ_total, Q_Term, Nr_Term
      integer :: icsf_nr
      integer :: I, N1_MAX, N1_big, N2_big, number_1, number_2
      integer :: I_Count, iterm_1, iterm_2, N1_big_jj, N2_big_jj
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      isubc_aviable = 3
      if(expansion_cfg_LS%csfs(1)%nosubc /= isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
           expansion_cfg_LS%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_cLSJ3 ', &
           'module transform_lsclsj3: cfg_Ji_structure%nr_of_subc.gt.',&
           'isubc_aviable'
         stop
      end if
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
      call define_number_of_csfs_cLSJ3(inr_of_csfs_cLSJ3)
      expansion_cfg_cLSJ3%size=inr_of_csfs_cLSJ3
      if(itype.ne.1) then
         allocate(expansion_cfg_cLSJ3%csfs(expansion_cfg_cLSJ3%size))
         allocate(expansion_cfg_cLSJ3%coeffs(expansion_cfg_cLSJ3%size))
         icsf_cLSJ3=0
         do icsf_LS=1, nomach_csfs_ls%list_size,1
           icsf_nr = nomach_csfs_ls%items(icsf_LS)
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il == 3 .and. &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big>=3) then
             Nr_Term = expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
             Q_Term = gettermLSQ(                                      &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,            &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,Nr_term,&
             expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,                &
             expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS)
           else
             Q_Term =                                                  &
                    2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1- &
                    expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
             Nr_Term = 0
           end if
           inosubc = expansion_cfg_LS%csfs(icsf_nr)%nosubc
           N1_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il
           if (N1_MAX == 0) then
              N1_big = 0
              N2_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              I_Count = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big   &
                                                        <= N1_MAX) then
              N1_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              N2_big= 0
              I_Count = N1_big + 1
           else
              N1_big=N1_MAX
              N2_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big-N1_big
              I_Count=                                                 &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+2      &
                -N2_big+1
           end if
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
           do L = L_min, L_max, 2
              do IS = IS_min, IS_max, 2
                 J_12_min = iabs(L - IS)
                 J_12_max = L + IS
                 do J_12 = J_12_min, J_12_max, 2
                    J_1_min = iabs(J_12 - iJ_total)
                    J_1_max = J_12 + iJ_total
                    do J_1 = J_1_min, J_1_max, 2
                       do I = 1,I_Count,1
                          if (N1_big-I+1 .ge. 0) then
                             if(icsf_cLSJ3.gt.expansion_cfg_cLSJ3%size) then
                               write(*,*)                           &
                              'stop at subroutine form_csfs_cLSJ3:',&
                               ' icsf_nr.gt.expansion_cfg_cLSJ3%size'
                               write(*,*)                           &
                               'J_1,icsf_cLSJ3', J_1,icsf_cLSJ3,    &
                               'expansion_cfg_cLSJ3%size',          &
                               expansion_cfg_cLSJ3%size
                               stop
                             end if
                            if(I == 1) then
                               if (N1_MAX /= 0) then
                                 N1_big_jj = N1_big
                                 call gettermjj                     &
                                 (2*expansion_cfg_LS%csfs(icsf_nr)% &
                                 subc_cfg(1)%il-1,N1_big_jj,        &
                                 jj_1_term,number_1)
                              else
                                 N1_big_jj = 0
                                 number_1 = 1
                                 jj_1_term(1)%j = 0
                                 jj_1_term(1)%Q = 0
                                 jj_1_term(1)%subshellJ = 0
                               end if
                               N2_big_jj = N2_big
                               call gettermjj                       &
                               (2*expansion_cfg_LS%csfs(icsf_nr)%   &
                               subc_cfg(1)%il+1,N2_big_jj,jj_2_term,&
                               number_2)
                             else
                               N1_big_jj = N1_big-I+1
                               N2_big_jj = N2_big+I-1
                               call gettermjj                       &
                               (2*expansion_cfg_LS%csfs(icsf_nr)%   &
                               subc_cfg(1)%il-1,N1_big_jj,jj_1_term,&
                               number_1)
                               call gettermjj                       &
                               (2*expansion_cfg_LS%csfs(icsf_nr)%   &
                               subc_cfg(1)%il+1,N2_big_jj,jj_2_term,&
                               number_2)
                             end if
                             do iterm_1 = 1,number_1
                             do iterm_2 = 1,number_2
                               if(ITTK(2*jj_1_term(iterm_1)%subshellJ, &
                               2*jj_2_term(iterm_2)%subshellJ,         &
                               2*J_1) == 0) cycle
!GG!                     coef = coefLSjj(                                  &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,    &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,&
!GG!                     Nr_Term,                                          &
!GG!                     Q_Term,                                           &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,        &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS,        &
!GG!                     J_1,                                              &
!GG!                     jj_1_term(iterm_1)%j,                             &
!GG!                     N1_big_jj,                                        &
!GG!                     jj_1_term(iterm_1)%Q,                             &
!GG!                     jj_1_term(iterm_1)%subshellJ,                     &
!GG!                     jj_2_term(iterm_2)%j,                             &
!GG!                     jj_2_term(iterm_2)%Q,                             &
!GG!                     jj_2_term(iterm_2)%subshellJ)
!GG!                                if (coef == ZERO_dp) cycle
                                icsf_cLSJ3 = icsf_cLSJ3 + 1
                                expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%nosubc= &
                                                                inosubc
                                allocate(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                subc_cfg(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                nosubc))
                                allocate(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                subc(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                nosubc))
                                allocate(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                iM1(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)% &
                                                               nosubc))
                                allocate(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                iM2(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)% &
                                                                nosubc))
                                allocate(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                iJ(expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%nosubc))
                                do isubc=1,expansion_cfg_cLSJ3%csfs(       &
                                                    icsf_cLSJ3)%nosubc,1
                                   expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%  &
                                   subc_cfg(isubc) = expansion_cfg_LS%csfs(&
                                   icsf_nr)%subc_cfg(isubc)
                                   expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%  &
                                   subc(isubc) = expansion_cfg_LS%csfs(    &
                                   icsf_nr)%subc(isubc)
                                   if(isubc == 1) then
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                     iM1(1) = N1_big_jj
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                     iM2(1) = N2_big_jj
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                                         iJ(1)=J_1
                                   else if(isubc == 2) then
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                                         iM1(2)=L
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                                        iM2(2)=IS
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                                        iJ(2)=J_12
                                  else if(isubc == 3) then
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                     iM1(3) = iterm_1
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                     iM2(3) = iterm_2
                                     expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%&
                                                          iJ(3)=iJ_total
                                  end if
                                end do
                             end do
                          end do
                          end if
                       end do
                    end do
                 end do
              end do
           end do
         end do
!
      end if
!
!write(*,*) '      end subroutine form_csfs_cLSJ3'
!write(iwrite_log,*) '      end subroutine form_csfs_cLSJ3'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine  define_number_of_csfs_cLSJ3(irez)
!                                                                      *
!     This subroutine defines the number of csfs in cLSJ3 coupling     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
         implicit none
         integer, intent(out) :: irez
!GG!         real(kind=dp) :: coef
         integer :: L, L_min, L_max, IS, IS_min, IS_max
         integer :: J_1, J_1_min, J_1_max
         integer :: J_12, J_12_min, J_12_max
         integer :: isubc_aviable
         integer :: icsf, icsf_nr, ITTK
         integer :: I, N1_MAX, N1_big, N2_big, number_1, number_2
         integer :: I_Count, iterm_1, iterm_2, N1_big_jj, N2_big_jj
         integer :: Q_Term, Nr_Term
         type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
         isubc_aviable = 3
         if(expansion_cfg_LS%csfs(1)%nosubc /= isubc_aviable) then
           write(*,*)'STOP at subroutine define_number_of_csfs_cLSJ3 ',&
            'module transform_lsclsJ3:cfg_Ji_structure%nr_of_subc.gt.',&
            'isubc_aviable'
           stop
         end if
         irez=0
         do icsf=1, nomach_csfs_ls%list_size,1
           icsf_nr = nomach_csfs_ls%items(icsf)
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il == 3 .and. &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big>=3) then
             Nr_Term = expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
             Q_Term = gettermLSQ(                                      &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,            &
             expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,Nr_term,&
             expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,                &
             expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS)
           else
             Q_Term =                                                  &
                    2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1- &
                    expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
             Nr_Term = 0
           end if
           N1_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il
           if (N1_MAX == 0) then
              N1_big = 0
              N2_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              I_Count = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big   &
                                                        <= N1_MAX) then
              N1_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              N2_big= 0
              I_Count = N1_big + 1
           else
              N1_big=N1_MAX
              N2_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big-N1_big
              I_Count=                                                 &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+2      &
                -N2_big+1
           end if
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
           do L = L_min, L_max, 2
              do IS = IS_min, IS_max, 2
                 J_12_min = iabs(L - IS)
                 J_12_max = L + IS
                 do J_12 = J_12_min, J_12_max, 2
                    J_1_min = iabs(J_12 - iJ_total)
                    J_1_max = J_12 + iJ_total
                    do J_1 = J_1_min, J_1_max, 2
                       do I = 1,I_Count,1
                          if (N1_big-I+1 .ge. 0) then
                             if(I == 1) then
                                if (N1_MAX /= 0) then
                                  N1_big_jj = N1_big
                                  call gettermjj(                      &
                                  2*expansion_cfg_LS%csfs(icsf_nr)%    &
                                  subc_cfg(1)%il-1,N1_big_jj,jj_1_term,&
                                  number_1)
                                else
                                  N1_big_jj = 0
                                  number_1 = 1
                                  jj_1_term(1)%j = 0
                                  jj_1_term(1)%Q = 0
                                  jj_1_term(1)%subshellJ = 0
                                end if
                                N2_big_jj = N2_big
                                call gettermjj(                        &
                                2*expansion_cfg_LS%csfs(icsf_nr)%      &
                                subc_cfg(1)%il+1,N2_big_jj,jj_2_term,  &
                                number_2)
                             else
                                N1_big_jj = N1_big-I+1
                                N2_big_jj = N2_big+I-1
                                call gettermjj(                        &
                                2*expansion_cfg_LS%csfs(icsf_nr)%      &
                                subc_cfg(1)%il-1,N1_big_jj,jj_1_term,  &
                                number_1)
                                call gettermjj(                        &
                                2*expansion_cfg_LS%csfs(icsf_nr)%      &
                                subc_cfg(1)%il+1,N2_big_jj,jj_2_term,  &
                                number_2)
                             end if
                             do iterm_1 = 1,number_1
                               do iterm_2 = 1,number_2
                                 if(ITTK(                              &
                                   2*jj_1_term(iterm_1)%subshellJ,     &
                                   2*jj_2_term(iterm_2)%subshellJ,     &
                                   2*J_1) == 0) cycle
!GG!                     coef = coefLSjj(                                  &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,    &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,&
!GG!                     Nr_Term,                                          &
!GG!                     Q_Term,                                           &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,        &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS,        &
!GG!                     J_1,                                              &
!GG!                     jj_1_term(iterm_1)%j,                             &
!GG!                     N1_big_jj,                                        &
!GG!                     jj_1_term(iterm_1)%Q,                             &
!GG!                     jj_1_term(iterm_1)%subshellJ,                     &
!GG!                     jj_2_term(iterm_2)%j,                             &
!GG!                     jj_2_term(iterm_2)%Q,                             &
!GG!                     jj_2_term(iterm_2)%subshellJ)
!GG!                                 if (coef /= ZERO_dp) irez = irez + 1
                                 irez = irez + 1
                               end do
                             end do
                          end if
                       end do
                    end do
                 end do
              end do
           end do
         end do
         end subroutine  define_number_of_csfs_cLSJ3
      end subroutine form_csfs_cLSJ3
!
!***********************************************************************
!                                                                      *
      subroutine matrix_LS_cLSJ3(icsf_LS,icsf_cLSJ3,rez)
!                                                                      *
!     This subroutine calculates the transformation                    *
!     matrix element between two csfs  (in LS and cLSJ3 couplings)     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
      integer,intent(in) :: icsf_LS, icsf_cLSJ3
      real(kind=dp),intent(out) :: rez
      real(kind=dp) :: rez_1, rez_2
      integer :: L_1, L_12, L_2, L_3, L_23, L
      integer :: IS_1, IS_12, IS_2, IS_3, IS_23, IS
      integer :: J_1, J_12, J
      integer :: nosubc, N1, N2, number_1, number_2, num_1, num_2
      integer :: Q_Term, Nr_Term
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      nosubc=expansion_cfg_LS%csfs(1)%nosubc
      rez = ZERO_dp
      if(equivalent_LScLSJ3(expansion_cfg_LS%csfs(icsf_LS),            &
         expansion_cfg_cLSJ3%csfs(icsf_cLSJ3))) then
         L_1    = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL
         L_12   = expansion_cfg_LS%csfs(icsf_LS)%iM1(2)
         L_2    = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iL
         L_3    = expansion_cfg_LS%csfs(icsf_LS)%subc(3)%iL
         L_23   = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM1(2)
         L      = expansion_cfg_LS%csfs(icsf_LS)%iM1(3)
         IS_1   = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS
         IS_12  = expansion_cfg_LS%csfs(icsf_LS)%iM2(2)
         IS_2   = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iS
         IS_3   = expansion_cfg_LS%csfs(icsf_LS)%subc(3)%iS
         IS_23  = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM2(2)
         IS     = expansion_cfg_LS%csfs(icsf_LS)%iM2(3)
         J_1    = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iJ(1)
         J_12   = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iJ(2)
         J      = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iJ(3)
         N1     = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM1(1)
         N2     = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM2(1)
         number_1 = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM1(3)
         number_2 = expansion_cfg_cLSJ3%csfs(icsf_cLSJ3)%iM2(3)
         if(expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il == 3 .and.   &
             expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big>=3) then
             Nr_Term = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%inr
             Q_Term = gettermLSQ(                                      &
             expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il,            &
             expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big,Nr_term,&
             expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL,                &
             expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS)
         else
             Q_Term =                                                  &
                    2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il+1- &
                    expansion_cfg_LS%csfs(icsf_LS)%subc(1)%inr
             Nr_Term = 0
         end if
         if (expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il /= 0) then
            call gettermjj                                             &
               (2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il-1,     &
               N1,jj_1_term, num_1)
         else
            jj_1_term(1)%j = 0
            jj_1_term(1)%Q = 0
            jj_1_term(1)%subshellJ = 0
         end if
         call gettermjj                                                &
            (2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il+1,        &
            N2,jj_2_term,num_2)
         rez_1 = coefLSjj(                                         &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il,    &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big,&
                 Nr_Term,                                          &
                 Q_Term,                                           &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL,        &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS,        &
                 J_1,                                              &
                 jj_1_term(number_1)%j,                            &
                 N1,                                               &
                 jj_1_term(number_1)%Q,                            &
                 jj_1_term(number_1)%subshellJ,                    &
                 jj_2_term(number_2)%j,                            &
                 jj_2_term(number_2)%Q,                            &
                 jj_2_term(number_2)%subshellJ)
         call cLSJ3PER(L_1,L_12, L_2, L_3, L_23, L,                 &
                       IS_1,IS_12,IS_2,IS_3,IS_23,IS,J_1,J_12,J,rez_2)
            rez = rez_1 * rez_2
      end if
!      print*, "isena is matrix_LS_cLSJ3"
      end subroutine matrix_LS_cLSJ3
!
!***********************************************************************
!                                                                      *
      subroutine dotransform_LScLSJ3
!                                                                      *
!     This subroutine calculates the weights of the                    *
!     expansions in the cLSJ3 coupling                                 *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp)::coeff_cLSJ3,coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_cLSJ3
!
      melement=0.
      sum_of_state = 0.
      do icsf_cLSJ3=1,expansion_cfg_cLSJ3%size
         coeff_cLSJ3 = 0.
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_cLSJ3(icsf_LS,icsf_cLSJ3,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_cLSJ3 = coeff_cLSJ3 + coeff_LS*melement
         end do
         if((dabs(coeff_cLSJ3)-dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_LScLSJ3:',&
               ' coeff_cLSJ3=',coeff_cLSJ3
            write(iwrite_log,*)'possible error at subroutine ',        &
               'dotransform_LScLSJ3: coeff_cLSJ3=',coeff_cLSJ3
         end if
         expansion_cfg_cLSJ3%coeffs(icsf_cLSJ3)= coeff_cLSJ3
         sum_of_state = sum_of_state + coeff_cLSJ3*coeff_cLSJ3
      end do
!
      if((sum_of_state-10*dp_coeff_precision).gt.ONE_dp) then
        write(*,*)'possible error at subroutine dotransform_LScLSJ3: ',&
            'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
            'dotransform_LScLSJ3: sum_of_state=',sum_of_state
      end if
!
      end subroutine dotransform_LScLSJ3
!
!***********************************************************************
!                                                                      *
      subroutine print_cfg_LScLSJ3(itype)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the unit "iwrite_cfg_expansions_cLSJ3"             *
!     (see module "Coupling_constants")                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in) :: itype
      integer             :: icsf, isubc, j
      integer             :: num_1, num_2, number_1, number_2
      character(len=1)    :: CVAL
      character(len=2)    :: LSJVAL
      character(len=4)    :: JVAL
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
!
!     LScLSJ3 - Coupling
!
      if(itype.gt.0) then
         write(iwrite_cfg_expansions_cLSJ3,*)  &
         '-------------------------------------'
         write(iwrite_cfg_expansions_cLSJ3,*)  &
         'state Nr.',expansion_cfg_LS%nr_of_state
         write(iwrite_cfg_expansions_cLSJ3,'(3x,a3,a4,a9,3x,f15.8)')   &
         'J =', JVAL(states%states(expansion_cfg_LS%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LS%nr_of_state)%energy
         write(iwrite_cfg_expansions_cLSJ3,'(3x,a30,i2)')              &
         'expansion size (LS coupling): ', expansion_cfg_LS%size
         write(iwrite_cfg_expansions_cLSJ3,*)''
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1,expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_cLSJ3,*)     &
                     'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_cLSJ3,                &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf,  &
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%in,       &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%iN_big,   &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_cLSJ3,                &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_cLSJ3,                &
                     '(1x,i5,3x,3(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,3),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
              'To many coupled shells in Coupling_transform_cfg_LScLSJ3'
                  end if
               else
                  write(iwrite_cfg_expansions_cLSJ3,'(9x,a51)')   &
                  'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LS%csfs(icsf)%subc)  &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM1) &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_cLSJ3,             &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_LS%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LS%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_cLSJ3,              &
                     '(12x,7(1x,i2,a1,i1))')                         &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,      &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL),&
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,       &
                      j=1,2),                                        &
                      expansion_cfg_LS%csfs(icsf)%iM2(2)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(2))
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_cLSJ3,               &
                     '(12x,3(1x,i2,a1,i1),3(1x,i2,a1))')              &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,       &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL), &
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,        &
                      j=1,3),                                         &
                      (expansion_cfg_LS%csfs(icsf)%iM2(j)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(j)),j=2,3)
                  end if
               else
                   write(iwrite_cfg_expansions_cLSJ3,'(9x,a63)') &
       'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_cLSJ3,*)' '
            end do
         else
            write(iwrite_cfg_expansions_cLSJ3,*) &
            'expansion_cfg_LS%csfs NOT associated'
         end if
      end if
!
!     cLSJ3 - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_cLSJ3,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_cLSJ3,*) &
         'state Nr.',expansion_cfg_cLSJ3%nr_of_state
         write(iwrite_cfg_expansions_cLSJ3,'(3x,a3,a4,a9,3x,f15.8)')   &
         'J =', JVAL(states%states(expansion_cfg_cLSJ3%nr_of_state)%J),&
       ' Energy =',states%states(expansion_cfg_cLSJ3%nr_of_state)%energy
         write(iwrite_cfg_expansions_cLSJ3,'(3x,a33,i2)')              &
         'expansion size (cLSJ3 coupling): ', expansion_cfg_cLSJ3%size
         write(iwrite_cfg_expansions_cLSJ3,*)''
         if(associated(expansion_cfg_cLSJ3%csfs)) then
            do icsf=1,expansion_cfg_cLSJ3%size
               if(associated(expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg)) then
!                  if(dabs(expansion_cfg_cLSJ3%coeffs(icsf)) <=          &
!                                              dp_coeff_precision) cycle
                  if (expansion_cfg_cLSJ3%csfs(icsf)%nosubc == 3) then
                     number_1 = expansion_cfg_cLSJ3%csfs(icsf)%iM1(3)
                     number_2 = expansion_cfg_cLSJ3%csfs(icsf)%iM2(3)
                     if(expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%il/=0)&
                                                                   then
                        call gettermjj                                 &
                          (2*expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%il&
                          -1,                                          &
                           expansion_cfg_cLSJ3%csfs(icsf)%iM1(1),      &
                           jj_1_term,num_1)
                     else
                        jj_1_term(number_1)%j = -1
                        jj_1_term(number_1)%subshellJ = 0
                        jj_1_term(number_1)%nu  = 0
                     end if
                     call gettermjj                                    &
                       (2*expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%il+1,&
                        expansion_cfg_cLSJ3%csfs(icsf)%iM2(1),         &
                        jj_2_term,num_2)
                     write(iwrite_cfg_expansions_cLSJ3,                &
                     '(1x,i5,3x,2(i2,a2,"(",i1,")"),2(i2,a1,"(",i1, &
                     ")"),6x,a6,f10.7)')icsf,                          &
                     expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%in,    &
                     LSJVAL(                                           &
                     expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%il,    &
                     jj_1_term(number_1)%j),                           &
                     expansion_cfg_cLSJ3%csfs(icsf)%iM1(1),            &
                     expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%in,    &
                     LSJVAL(                                           &
                     expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(1)%il,    &
                     jj_2_term(number_2)%j),                           &
                     expansion_cfg_cLSJ3%csfs(icsf)%iM2(1),            &
                     (expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(j)%in,   &
                 CVAL(1,expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg(j)%iN_big,&
                     j=2,3),                                           &
                     'coeff:',expansion_cfg_cLSJ3%coeffs(icsf)
                  end if
               else
                  write(iwrite_cfg_expansions_cLSJ3,'(9x,a51)')        &
                'expansion_cfg_cLSJ3%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_cLSJ3%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_cLSJ3%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_cLSJ3%csfs(icsf)%iM2)) then
                  if (expansion_cfg_cLSJ3%csfs(icsf)%nosubc == 3) then
                     write(iwrite_cfg_expansions_cLSJ3,                &
                     '(10x,2(1x,a4,i2)," (",i2,a1,i1,1x,i2,a1,i1,")",  &
                     i2,a1,1x," [",a4,","a4," ]")')                   &
                     JVAL(jj_1_term(number_1)%subshellJ),              &
                     jj_1_term(number_1)%nu,                           &
                     JVAL(jj_2_term(number_2)%subshellJ),              &
                     jj_2_term(number_2)%nu,                           &
                     (expansion_cfg_cLSJ3%csfs(icsf)%subc(j)%iS+1,     &
                     CVAL(2,expansion_cfg_cLSJ3%csfs(icsf)%subc(j)%iL),&
                     expansion_cfg_cLSJ3%csfs(icsf)%subc(j)%inr,       &
                     j=2,3),                                           &
                     expansion_cfg_cLSJ3%csfs(icsf)%iM2(2)+1,          &
                     CVAL(2,expansion_cfg_cLSJ3%csfs(icsf)%iM1(2)),    &
                     JVAL(expansion_cfg_cLSJ3%csfs(icsf)%iJ(1)),       &
                     JVAL(expansion_cfg_cLSJ3%csfs(icsf)%iJ(2))
                  end if
               else
                  write(iwrite_cfg_expansions_cLSJ3,'(9x,a63)')        &
                     'expansion_cfg_cLSJ3%csfs(icsf)%subc (',          &
                     'or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_cLSJ3,*)' '
            end do
         else
            write(iwrite_cfg_expansions_cLSJ3,*) &
            'expansion_cfg_cLSJ3%csfs NOT associated'
         end if
            write(iwrite_cfg_expansions_cLSJ3,*) &
            '-------------------------------------'
      end if
      end subroutine print_cfg_LScLSJ3
!
      end module Coupling_transform_cfg_LScLSJ3
