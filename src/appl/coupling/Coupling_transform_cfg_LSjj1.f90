!
!***********************************************************************
!
      module Coupling_transform_cfg_LSjj1
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
      public  :: main_cfg_lslsjj1
      public  :: count_nr_of_csfs_jj1
      public  :: delete_cfg_expansions
      private :: form_list_nomach_csfs
      private :: form_csfs_jj1
      private :: matrix_LS_jj1
      private :: dotransform_LSjj1
      private :: print_cfg_LSjj1
      private :: equivalent_LSjj1
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
      type(expansion), public :: expansion_cfg_jj1
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_cfg_lslsjj1(print_level)
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
!      write(*,*) '      subroutine main_cfg_lslsjj1'
      expansion_cfg_jj1%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=0
      call form_csfs_jj1(itype)
      call dotransform_LSjj1
      if(print_level.gt.1) call print_cfg_LSjj1(2) ! 2-means print LS and jj1 expansions
      deallocate(nomach_csfs_ls%items, STAT = error)
!      write(*,*) '      end subroutine main_cfg_lslsjj1'
      end subroutine main_cfg_lslsjj1
!
!***********************************************************************
!                                                                      *
      subroutine count_nr_of_csfs_jj1(nr_of_csfs)
!                                                                      *
!     This subroutine counts the number of oneconfigurational          *
!     expansion in jj1 coupling                                        *
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
      if(expansion_cfg_LS%csfs(1)%nosubc /= 1) then
         write(*,*)                                                    &
            'ERROR at subroutine count_nr_of_csfs_jj1: nosubc.lt.2'
         write(*,*)                                                    &
           '(you can not use jj1 coupling with bigger than 1 subshells)'
         write(*,*)'program is terminated'
         stop
      end if
      expansion_cfg_jj1%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=1
      call form_csfs_jj1(itype)
      nr_of_csfs = expansion_cfg_jj1%size
      if(associated(nomach_csfs_ls%items))                             &
                             deallocate(nomach_csfs_ls%items,STAT=error)
      end subroutine count_nr_of_csfs_jj1
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_expansions
!                                                                      *
!     This subroutine deallocates the arrays of                        *
!     oneconfigurational expansions "expansion_cfg_jj1" and            *
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
!delete expansion_cfg_jj1
      if(associated(expansion_cfg_jj1%coeffs))                         &
                                      deallocate(expansion_cfg_jj1%csfs)
      if(associated(expansion_cfg_jj1%csfs)) then
         do icsf=1, expansion_cfg_jj1%size
            if(associated(expansion_cfg_jj1%csfs(icsf)%subc_cfg))      &
                       deallocate(expansion_cfg_jj1%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_jj1%csfs(icsf)%subc))          &
                           deallocate(expansion_cfg_jj1%csfs(icsf)%subc)
            if(associated(expansion_cfg_jj1%csfs(icsf)%iM1))           &
                            deallocate(expansion_cfg_jj1%csfs(icsf)%iM1)
            if(associated(expansion_cfg_jj1%csfs(icsf)%iM2))           &
                            deallocate(expansion_cfg_jj1%csfs(icsf)%iM2)
            if(associated(expansion_cfg_jj1%csfs(icsf)%iJ))            &
                             deallocate(expansion_cfg_jj1%csfs(icsf)%iJ)
        end do
        deallocate(expansion_cfg_jj1%csfs)
      end if
      expansion_cfg_jj1%size=0
      expansion_cfg_jj1%nr_of_state=0
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
      function equivalent_LSjj1(csf1, csf2)   result(rez)
!                                                                      *
!     This function defines the "equivalency" of                       *
!     two csfs in LS coupling for the formation of jj1 csfs            *
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
      if(csf1%nosubc /= 1 .or.csf2%nosubc /= 1) then
        write(*,*)'ERROR at subroutine equivalent_LSjj1: ',            &
                   'csf1%nosubc.lt.2.or.csf2%nosubc.lt.2'
        write(*,*)'(you can not use jj1 coupling with bigger than 1 ', &
                   'subshells)'
         write(*,*)'program will be terminated'
         stop
      end if
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
!            else if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
!               if(isubc /= 2) rez=.false.
!               rez=.false.
!               exit
            end if
         end do
      else
         rez=.false.
      end if
      end function equivalent_LSjj1
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
            if(equivalent_LSjj1(expansion_cfg_LS%csfs(icsf_LS),        &
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
      subroutine form_csfs_jj1(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in jj1 coupling "expansion_cfg_jj1",                   *
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
                      !of csfs_jj1 only,
                      !itype.ne.1 - perform full
                      !calculation
!GG!      real(kind=dp) :: coef
      integer :: inr_of_csfs_jj1
      integer :: icsf_LS, icsf_jj1
      integer :: inosubc, isubc, isubc_aviable
      integer :: icsf, icsf_nr
      integer :: I, N1_MAX, ITTK
      integer :: N1_big, N2_big
      integer :: number_1, number_2
      integer :: I_Count1
      integer :: iterm_1, iterm_2
      integer :: N1_big_jj, N2_big_jj
      integer :: iJ_total, Q_Term1, Nr_Term1
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
!
      isubc_aviable = 1
      if(expansion_cfg_LS%csfs(1)%nosubc.gt.isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_jj1 ',   &
            'module transform_lsjj: cfg_Ji_structure%nr_of_subc.gt.',  &
            'isubc_aviable'
         stop
      end if
!
!      write(*,*) '      subroutine form_csfs_jj1'
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
!
      call define_number_of_csfs_jj1(inr_of_csfs_jj1)
!
      expansion_cfg_jj1%size=inr_of_csfs_jj1
!
      if(itype.ne.1) then
         allocate(expansion_cfg_jj1%csfs(expansion_cfg_jj1%size))
         allocate(expansion_cfg_jj1%coeffs(expansion_cfg_jj1%size))
         icsf_jj1=0
         do icsf_LS=1, nomach_csfs_ls%list_size,1
            icsf_nr = nomach_csfs_ls%items(icsf_LS)
!
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il == 3 .and. &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big>=3) then
            Nr_Term1 = expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
            Q_Term1 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,             &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,Nr_term1,&
            expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,                 &
            expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS)
           else
            Q_Term1 =                                                  &
                   2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1-  &
                   expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
            Nr_Term1 = 0
           end if
           inosubc = expansion_cfg_LS%csfs(icsf_nr)%nosubc
           N1_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il
           if (N1_MAX == 0) then
              N1_big = 0
              N2_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              I_Count1 = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big   &
                                                        <= N1_MAX) then
              N1_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              N2_big= 0
              I_Count1 = N1_big + 1
           else
              N1_big=N1_MAX
              N2_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big-N1_big
              I_Count1 =                                               &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+2      &
                -N2_big+1
           end if
!
           do I = 1,I_Count1,1
              if (N1_big-I+1 .ge. 0) then
                if(I == 1) then
                   if (N1_MAX /= 0) then
                      N1_big_jj = N1_big
                      call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)% &
                         subc_cfg(1)%il-1,N1_big_jj,jj_1_term,number_1)
                   else
                      N1_big_jj = 0
                      number_1 = 1
                      jj_1_term(1)%j = 0
                      jj_1_term(1)%Q = 0
                      jj_1_term(1)%subshellJ = 0
                   end if
                   N2_big_jj = N2_big
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il+1,N2_big_jj,jj_2_term,number_2)
                else
                   N1_big_jj = N1_big-I+1
                   N2_big_jj = N2_big+I-1
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il-1,N1_big_jj,jj_1_term,number_1)
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il+1,N2_big_jj,jj_2_term,number_2)
                end if
                do iterm_1 = 1,number_1
                  do iterm_2 = 1,number_2
                     if(ITTK(2*jj_1_term(iterm_1)%subshellJ,           &
                         2*jj_2_term(iterm_2)%subshellJ,               &
                                                2*iJ_total) == 0) cycle
!GG!                     coef = coefLSjj(                                  &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,    &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,&
!GG!                     Nr_Term1,                                         &
!GG!                     Q_Term1,                                          &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,        &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS,        &
!GG1                     iJ_total,                                         &
!GG!                     jj_1_term(iterm_1)%j,                             &
!GG!                     N1_big_jj,                                        &
!GG!                     jj_1_term(iterm_1)%Q,                             &
!GG!                     jj_1_term(iterm_1)%subshellJ,                     &
!GG!                     jj_2_term(iterm_2)%j,                             &
!GG!                     jj_2_term(iterm_2)%Q,                             &
!GG!                     jj_2_term(iterm_2)%subshellJ)
!GG!                     if (coef == ZERO_dp) cycle
!
                     icsf_jj1 = icsf_jj1 + 1
                     expansion_cfg_jj1%csfs(icsf_jj1)%nosubc=inosubc
                     allocate(expansion_cfg_jj1%csfs(icsf_jj1)%        &
                          subc_cfg(1))
                     allocate(expansion_cfg_jj1%csfs(icsf_jj1)%subc(1))
                     allocate(expansion_cfg_jj1%csfs(icsf_jj1)%iM1(1))
                     allocate(expansion_cfg_jj1%csfs(icsf_jj1)%iM2(1))
                     allocate(expansion_cfg_jj1%csfs(icsf_jj1)%iJ(1))
                     expansion_cfg_jj1%csfs(icsf_jj1)%iJ(1) = 0
                     expansion_cfg_jj1%csfs(icsf_jj1)%iM1(1)           &
                                           = 1000*N1_big_jj + N2_big_jj
                     expansion_cfg_jj1%csfs(icsf_jj1)%iM2(1)           &
                                           = 1000*iterm_1 + iterm_2
                     expansion_cfg_jj1%csfs(icsf_jj1)%subc_cfg(1)=     &
                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)
!
                     expansion_cfg_jj1%csfs(icsf_jj1)%subc(1) =        &
                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)
                  end do
                end do
              end if
           end do
         end do
!
      end if
!      write(*,*) '      end subroutine form_csfs_jj1'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine  define_number_of_csfs_jj1(irez)
!                                                                      *
!     This subroutine defines the number of csfs in jj1 coupling       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
         implicit none
         integer, intent(out)::irez
!GG!         real(kind=dp) :: coef
         integer :: inosubc, isubc, isubc_aviable
         integer :: icsf, icsf_nr, ITTK, ITREXG
         integer :: I, N1_MAX
         integer :: N1_big, N2_big
         integer :: number_1, number_2
         integer :: I_Count1
         integer :: iterm_1, iterm_2
         integer :: N1_big_jj, N2_big_jj
         integer :: Q_Term1, Nr_Term1, Q_Term2, Nr_Term2
         type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
!      write(*,*) '      subroutine define_number_of_csfs_jj1'
         isubc_aviable = 1
         if(expansion_cfg_LS%csfs(1)%nosubc .gt. isubc_aviable) then
           write(*,*) 'STOP at subroutine define_number_of_csfs_jj1 ', &
              'module transform_lsjj: cfg_Ji_structure%nr_of_subc.gt.',&
              'isubc_aviable'
           stop
         end if
         irez=0
         do icsf=1, nomach_csfs_ls%list_size,1
            icsf_nr = nomach_csfs_ls%items(icsf)
!
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il == 3 .and. &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big>=3) then
            Nr_Term1 = expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
            Q_Term1 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,             &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,Nr_term1,&
            expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,                 &
            expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS)
           else
            Q_Term1 =                                                  &
                    2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1- &
                    expansion_cfg_LS%csfs(icsf_nr)%subc(1)%inr
            Nr_Term1 = 0
           end if
           inosubc = expansion_cfg_LS%csfs(icsf_nr)%nosubc
           N1_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il
           if (N1_MAX == 0) then
              N1_big = 0
              N2_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              I_Count1 = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big   &
                                                        <= N1_MAX) then
              N1_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big
              N2_big= 0
              I_Count1 = N1_big + 1
           else
              N1_big=N1_MAX
              N2_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big-N1_big
              I_Count1 =                                               &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+2      &
                -N2_big+1
           end if
!
           do I = 1,I_Count1,1
              if (N1_big-I+1 .ge. 0) then
                if(I == 1) then
                   if (N1_MAX /= 0) then
                      N1_big_jj = N1_big
                      call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)% &
                         subc_cfg(1)%il-1,N1_big_jj,jj_1_term,number_1)
                   else
                      N1_big_jj = 0
                      number_1 = 1
                      jj_1_term(1)%j = 0
                      jj_1_term(1)%Q = 0
                      jj_1_term(1)%subshellJ = 0
                   end if
                   N2_big_jj = N2_big
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il+1,N2_big_jj,jj_2_term,number_2)
                else
                   N1_big_jj = N1_big-I+1
                   N2_big_jj = N2_big+I-1
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il-1,N1_big_jj,jj_1_term,number_1)
                   call gettermjj(2*expansion_cfg_LS%csfs(icsf_nr)%    &
                         subc_cfg(1)%il+1,N2_big_jj,jj_2_term,number_2)
                end if
                do iterm_1 = 1,number_1
                  do iterm_2 = 1,number_2
                     if(ITTK(2*jj_1_term(iterm_1)%subshellJ,           &
                         2*jj_2_term(iterm_2)%subshellJ,               &
                                                2*iJ_total) == 0) cycle
!GG!                     coef = coefLSjj(                                  &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,    &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,&
!GG!                     Nr_Term1,                                         &
!GG!                     Q_Term1,                                          &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iL,        &
!GG!                     expansion_cfg_LS%csfs(icsf_nr)%subc(1)%iS,        &
!GG!                     iJ_total,                                         &
!GG!                     jj_1_term(iterm_1)%j,                             &
!GG!                     N1_big_jj,                                        &
!GG!                     jj_1_term(iterm_1)%Q,                             &
!GG!                     jj_1_term(iterm_1)%subshellJ,                     &
!GG!                     jj_2_term(iterm_2)%j,                             &
!GG!                     jj_2_term(iterm_2)%Q,                             &
!GG!                     jj_2_term(iterm_2)%subshellJ)
!GG!                     if (coef == ZERO_dp) cycle
!
                     irez = irez +1
                  end do
                end do
              end if
           end do
         end do
!      write(*,*) '     end subroutine define_number_of_csfs_jj1'
         end subroutine  define_number_of_csfs_jj1
      end subroutine form_csfs_jj1
!
!***********************************************************************
!                                                                      *
      subroutine matrix_LS_jj1(icsf_LS,icsf_jj1,rez)
!                                                                      *
!     This subroutine calculates the transformation                    *
!     matrix element between two csfs  (in LS and jj1 couplings)       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
      integer,intent(in):: icsf_LS,icsf_jj1
      real(kind=dp),intent(out)::rez
      integer :: nosubc, N1, N2
      integer :: J, number_1, number_2
      integer :: num_1, num_2
      integer :: Q_Term1, Nr_Term1
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
!      write(*,*) '      subroutine matrix_LS_jj1'
      nosubc=expansion_cfg_LS%csfs(1)%nosubc
      rez = ZERO_dp
      if(equivalent_LSjj1(expansion_cfg_LS%csfs(icsf_LS),            &
         expansion_cfg_jj1%csfs(icsf_jj1))) then
         J   = states%states(expansion_cfg_LS%nr_of_state)%J
         N1  = JTHN(expansion_cfg_jj1%csfs(icsf_jj1)%iM1(1),2,1000)
         N2  = JTHN(expansion_cfg_jj1%csfs(icsf_jj1)%iM1(1),1,1000)
         number_1 = JTHN(expansion_cfg_jj1%csfs(icsf_jj1)%iM2(1),2,1000)
         number_2 = JTHN(expansion_cfg_jj1%csfs(icsf_jj1)%iM2(1),1,1000)
         if(expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il == 3 .and. &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big>=3) then
            Nr_Term1 = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%inr
            Q_Term1 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il,             &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big,Nr_term1,&
            expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL,                 &
            expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS)
         else
            Q_Term1 =                                                  &
                   2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il+1-  &
                   expansion_cfg_LS%csfs(icsf_LS)%subc(1)%inr
            Nr_Term1 = 0
         end if
!
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
!
         rez   = coefLSjj(                                         &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il,    &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big,&
                 Nr_Term1,                                         &
                 Q_Term1,                                          &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL,        &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS,        &
                 J,                                                &
                 jj_1_term(number_1)%j,                            &
                 N1,                                               &
                 jj_1_term(number_1)%Q,                            &
                 jj_1_term(number_1)%subshellJ,                    &
                 jj_2_term(number_2)%j,                            &
                 jj_2_term(number_2)%Q,                            &
                 jj_2_term(number_2)%subshellJ)
      end if
!      write(*,*) '      end subroutine matrix_LS_jj1'
      end subroutine matrix_LS_jj1
!
!***********************************************************************
!                                                                      *
      subroutine dotransform_LSjj1
!                                                                      *
!     This subroutine calculates the weights of the                    *
!     expansions in the jj1 coupling                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp)::coeff_jj1,coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_jj1
!      write(*,*) '      subroutine dotransform_LSjj1'
      melement=0.
      sum_of_state = 0.
      do icsf_jj1=1,expansion_cfg_jj1%size
         coeff_jj1 = 0.
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_jj1(icsf_LS,icsf_jj1,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_jj1 = coeff_jj1 + coeff_LS*melement
         end do
         if((dabs(coeff_jj1)-dp_coeff_precision).gt.ONE_dp) then
           write(*,*)'possible error at subroutine dotransform_LSjj1:',&
               ' coeff_jj1=',coeff_jj1
           write(iwrite_log,*)'possible error at subroutine ',         &
               'dotransform_LSjj1: coeff_jj1=',coeff_jj1
         end if
         expansion_cfg_jj1%coeffs(icsf_jj1)= coeff_jj1
         sum_of_state = sum_of_state + coeff_jj1*coeff_jj1
      end do
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_LSjj1: ', &
            'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
            'dotransform_LSjj1: sum_of_state=',sum_of_state
      end if
!      write(*,*) '      end subroutine dotransform_LSjj1'
      end subroutine dotransform_LSjj1
!
!***********************************************************************
!                                                                      *
      subroutine print_cfg_LSjj1(itype)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the unit "iwrite_cfg_expansions_jj1"               *
!     (see module "Coupling_constants")                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in) :: itype
      integer             :: icsf, isubc, j!, JTHN
      integer             :: num_1, num_2
      integer             :: number_1, number_2
      character(len=1)    :: CVAL
      character(len=2)    :: LSJVAL
      character(len=4)    :: JVAL
      type(subshell_term), dimension(1:63) :: jj_1_term, jj_2_term
!
!     LS - Coupling
!
!      write(*,*) '      subroutine print_cfg_LSjj1'
      if(itype.gt.0) then
         write(iwrite_cfg_expansions_jj123,*)  &
         '-------------------------------------'
         write(iwrite_cfg_expansions_jj123,*)  &
         'state Nr.',expansion_cfg_LS%nr_of_state
         write(iwrite_cfg_expansions_jj123,'(3x,a3,a4,a9,3x,f15.8)')   &
         'J =', JVAL(states%states(expansion_cfg_LS%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LS%nr_of_state)%energy
         write(iwrite_cfg_expansions_jj123,'(3x,a29,i2)')              &
         'expansion size (LS coupling): ', expansion_cfg_LS%size
         write(iwrite_cfg_expansions_jj123,*)''
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1,expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
                  if(icsf==1)write(iwrite_cfg_expansions_jj123,*)      &
                                                              'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_jj123,                &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf,  &
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%in,       &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%iN_big,   &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
                 'To many couled shells in Coupling_transform_cfg_LSjj1'
                  end if
               else
                  write(iwrite_cfg_expansions_jj123,'(9x,a51)')   &
                  'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LS%csfs(icsf)%subc)  &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM1) &
               .and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_jj123,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_LS%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LS%csfs(icsf)%subc(1)%inr
                  end if
               else
                   write(iwrite_cfg_expansions_jj123,'(9x,a63)') &
       'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_jj123,*)' '
            end do
         else
            write(iwrite_cfg_expansions_jj123,*) &
            'expansion_cfg_LS%csfs NOT associated'
         end if
      end if
!
!     jj1 - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_jj123,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_jj123,*) &
         'state Nr.',expansion_cfg_jj1%nr_of_state
         write(iwrite_cfg_expansions_jj123,'(3x,a3,a4,a9,3x,f15.8)')   &
         'J =', JVAL(states%states(expansion_cfg_jj1%nr_of_state)%J),  &
         ' Energy =',states%states(expansion_cfg_jj1%nr_of_state)%energy
         write(iwrite_cfg_expansions_jj123,'(3x,a29,i2)')              &
         'expansion size (jj1 coupling): ', expansion_cfg_jj1%size
         write(iwrite_cfg_expansions_jj123,*)''
         if(associated(expansion_cfg_jj1%csfs)) then
            do icsf=1,expansion_cfg_jj1%size
               if(associated(expansion_cfg_jj1%csfs(icsf)%subc_cfg)) then
                  if (expansion_cfg_jj1%csfs(icsf)%nosubc == 1) then
                     number_1 =                                        &
                       JTHN(expansion_cfg_jj1%csfs(icsf)%iM2(1),2,1000)
                     number_2 =                                        &
                       JTHN(expansion_cfg_jj1%csfs(icsf)%iM2(1),1,1000)
                     if(expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%il/=0)&
                                                                   then
                        call gettermjj                                 &
                        (2*expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%il &
                        -1,                                            &
                      JTHN(expansion_cfg_jj1%csfs(icsf)%iM1(1),2,1000),&
                        jj_1_term,num_1)
                     else
                        jj_1_term(number_1)%j = -1
                        jj_1_term(number_1)%subshellJ = 0
                        jj_1_term(number_1)%nu  = 0
                     end if
                     call gettermjj                                    &
                     (2*expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%il+1, &
                     JTHN(expansion_cfg_jj1%csfs(icsf)%iM1(1),1,1000), &
                     jj_2_term,num_2)
                     write(iwrite_cfg_expansions_jj123,                &
                     '(1x,i5,3x,i2,a2,"(",i1,")",                      &
                     i2,a2,"(",i1,")",6x,a6,f10.7)')                   &
                     icsf,                                             &
                     expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%il,      &
                     jj_1_term(number_1)%j),                           &
                     JTHN(expansion_cfg_jj1%csfs(icsf)%iM1(1),2,1000), &
                     expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj1%csfs(icsf)%subc_cfg(1)%il,      &
                     jj_2_term(number_2)%j),                           &
                     JTHN(expansion_cfg_jj1%csfs(icsf)%iM1(1),1,1000), &
                     'coeff:',expansion_cfg_jj1%coeffs(icsf)
                  end if
               else
                  write(iwrite_cfg_expansions_jj123,'(9x,a51)') &
                  'expansion_cfg_jj1%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_jj1%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_jj1%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_jj1%csfs(icsf)%iM2)) then
                  if (expansion_cfg_jj1%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_jj123,                &
                     '(9x,2(a4,1x,i1,1x),"[",a4,"]")')                 &
                     JVAL(jj_1_term(number_1)%subshellJ),              &
                     jj_1_term(number_1)%nu,                           &
                     JVAL(jj_2_term(number_2)%subshellJ),              &
                     jj_2_term(number_2)%nu,                           &
                    JVAL(states%states(expansion_cfg_jj1%nr_of_state)%J)
                  end if
               else
                  write(iwrite_cfg_expansions_jj123,'(9x,a63)') &
       'expansion_cfg_jj1%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_jj123,*)' '
            end do
         else
            write(iwrite_cfg_expansions_jj123,*) &
            'expansion_cfg_jj1%csfs NOT associated'
         end if
            write(iwrite_cfg_expansions_jj123,*) &
            '-------------------------------------'
      end if
!      write(*,*) '      end subroutine print_cfg_LSjj1'
      end subroutine print_cfg_LSjj1
!
      end module Coupling_transform_cfg_LSjj1
