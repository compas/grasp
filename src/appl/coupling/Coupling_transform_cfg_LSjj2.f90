!
!***********************************************************************
!
      module Coupling_transform_cfg_LSjj2
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
      public  :: main_cfg_lslsjj2
      public  :: count_nr_of_csfs_jj2
      public  :: delete_cfg_expansions
      private :: form_list_nomach_csfs
      private :: form_csfs_jj2
      private :: matrix_LS_jj2
      private :: dotransform_LSjj2
      private :: print_cfg_LSjj2
      private :: equivalent_LSjj2
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
      type(expansion), public :: expansion_cfg_jj2
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_cfg_lslsjj2(print_level)
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
!      write(*,*) '      subroutine main_cfg_lslsjj2'
      expansion_cfg_jj2%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=0
      call form_csfs_jj2(itype)
      call dotransform_LSjj2
      if(print_level.gt.1) call print_cfg_LSjj2(2) ! 2-means print LS and jj2 expansions
      deallocate(nomach_csfs_ls%items, STAT = error)
!      write(*,*) '      end subroutine main_cfg_lslsjj2'
      end subroutine main_cfg_lslsjj2
!
!***********************************************************************
!                                                                      *
      subroutine count_nr_of_csfs_jj2(nr_of_csfs)
!                                                                      *
!     This subroutine counts the number of oneconfigurational          *
!     expansion in jj2 coupling                                        *
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
            'ERROR at subroutine count_nr_of_csfs_jj2: nosubc.lt.2'
         write(*,*)                                                    &
            '(you can not use jj2 coupling with less than 2 subshells)'
         write(*,*)'program is terminated'
         stop
      end if
      expansion_cfg_jj2%nr_of_state=expansion_cfg_LS%nr_of_state
      call form_list_nomach_csfs
      itype=1
      call form_csfs_jj2(itype)
      nr_of_csfs = expansion_cfg_jj2%size
      if(associated(nomach_csfs_ls%items))                             &
                             deallocate(nomach_csfs_ls%items,STAT=error)
      end subroutine count_nr_of_csfs_jj2
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_expansions
!                                                                      *
!     This subroutine deallocates the arrays of                        *
!     oneconfigurational expansions "expansion_cfg_jj2" and            *
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
!delete expansion_cfg_jj2
      if(associated(expansion_cfg_jj2%coeffs))                         &
                                      deallocate(expansion_cfg_jj2%csfs)
      if(associated(expansion_cfg_jj2%csfs)) then
         do icsf=1, expansion_cfg_jj2%size
            if(associated(expansion_cfg_jj2%csfs(icsf)%subc_cfg))      &
                       deallocate(expansion_cfg_jj2%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_jj2%csfs(icsf)%subc))          &
                           deallocate(expansion_cfg_jj2%csfs(icsf)%subc)
            if(associated(expansion_cfg_jj2%csfs(icsf)%iM1))           &
                            deallocate(expansion_cfg_jj2%csfs(icsf)%iM1)
            if(associated(expansion_cfg_jj2%csfs(icsf)%iM2))           &
                            deallocate(expansion_cfg_jj2%csfs(icsf)%iM2)
            if(associated(expansion_cfg_jj2%csfs(icsf)%iJ))            &
                             deallocate(expansion_cfg_jj2%csfs(icsf)%iJ)
        end do
        deallocate(expansion_cfg_jj2%csfs)
      end if
      expansion_cfg_jj2%size=0
      expansion_cfg_jj2%nr_of_state=0
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
      function equivalent_LSjj2(csf1, csf2)   result(rez)
!                                                                      *
!     This function defines the "equivalency" of                       *
!     two csfs in LS coupling for the formation of jj2 csfs             *
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
        write(*,*)'ERROR at subroutine equivalent_LSjj2: ',            &
                   'csf1%nosubc.lt.2.or.csf2%nosubc.lt.2'
        write(*,*)'(you can not use jj2 coupling with less than 2 ',   &
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
      end function equivalent_LSjj2
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
            if(equivalent_LSjj2(expansion_cfg_LS%csfs(icsf_LS),        &
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
      subroutine form_csfs_jj2(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in jj2 coupling "expansion_cfg_jj2",                   *
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
                      !of csfs_jj2 only,
                      !itype.ne.1 - perform full
                      !calculation
!GG!      real(kind=dp) :: coef
      integer :: inr_of_csfs_jj2
      integer :: icsf_LS, icsf_jj2
      integer :: J_1, J_1_min, J_1_max
      integer :: J_2, J_2_min, J_2_max
      integer :: IKK,JJTM,JJTP,JJTM_min,JJTM_max,JJTP_min,JJTP_max
      integer :: inosubc, isubc, isubc_aviable
      integer :: icsf, icsf_nr, ITTK, ITREXG
      integer :: I, II, N1_MAX, N2_MAX
      integer :: N1_big, N2_big, N3_big, N4_big
      integer :: number_1, number_2, number_3, number_4
      integer :: I_Count1, I_Count2
      integer :: iterm_1, iterm_2, iterm_3, iterm_4
      integer :: N1_big_jj, N2_big_jj, N3_big_jj, N4_big_jj
      integer :: iJ_total, Q_Term1, Nr_Term1, Q_Term2, Nr_Term2
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
!
      isubc_aviable = 2
      if(expansion_cfg_LS%csfs(1)%nosubc.gt.isubc_aviable) then
         write(*,*) "is available =",isubc_aviable," we have =",       &
            expansion_cfg_LS%csfs(1)%nosubc
         write(*,*) 'STOP at subroutine define_number_of_csfs_jj2 ',   &
            'module transform_lsjj: cfg_Ji_structure%nr_of_subc.gt.',  &
            'isubc_aviable'
         stop
      end if
!
!      write(*,*) '      subroutine form_csfs_jj2'
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
!
      call define_number_of_csfs_jj2(inr_of_csfs_jj2)
!
      expansion_cfg_jj2%size=inr_of_csfs_jj2
!
      if(itype.ne.1) then
         allocate(expansion_cfg_jj2%csfs(expansion_cfg_jj2%size))
         allocate(expansion_cfg_jj2%coeffs(expansion_cfg_jj2%size))
         icsf_jj2=0
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
              Q_Term1 =                                                &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1-&
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
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il == 3 .and. &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big>=3) then
            Nr_Term2 = expansion_cfg_LS%csfs(icsf_nr)%subc(2)%inr
            Q_Term2 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il,             &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big,Nr_term2,&
            expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL,                 &
            expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS)
           else
            Q_Term2 =                                                  &
                   2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+1-  &
                   expansion_cfg_LS%csfs(icsf_nr)%subc(2)%inr
            Nr_Term2 = 0
           end if
           N2_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il
           if (N2_MAX == 0) then
              N3_big = 0
              N4_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big
              I_Count2 = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big   &
                                                        <= N2_MAX) then
              N3_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big
              N4_big = 0
              I_Count2 = N3_big + 1
           else
              N3_big=N2_MAX
              N4_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big-N3_big
              I_Count2 =                                               &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+2      &
                -N4_big+1
           end if
!
           call getJMinMax(                                            &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,         &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,     &
                J_1_min,J_1_max)
           do J_1 = J_1_min, J_1_max, 2
                 do I = 1,I_Count1,1
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
!GG!                     Nr_Term1,                                         &
!GG!                     Q_Term1,                                          &
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
!GG!                          if (coef == ZERO_dp) cycle
!
                 do II = 1,I_Count2,1
                    if (N3_big-II+1 .ge. 0) then
                      if(II == 1) then
                        if (N2_MAX /= 0) then
                           N3_big_jj = N3_big
                           call gettermjj(                      &
                           2*expansion_cfg_LS%csfs(icsf_nr)%    &
                           subc_cfg(2)%il-1,N3_big_jj,jj_3_term,&
                           number_3)
                        else
                           N3_big_jj = 0
                           number_3 = 1
                           jj_3_term(1)%j = 0
                           jj_3_term(1)%Q = 0
                           jj_3_term(1)%subshellJ = 0
                        end if
                        N4_big_jj = N4_big
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il+1,N4_big_jj,jj_4_term,  &
                        number_4)
                      else
                        N3_big_jj = N3_big-II+1
                        N4_big_jj = N4_big+II-1
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il-1,N3_big_jj,jj_3_term,  &
                        number_3)
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il+1,N4_big_jj,jj_4_term,  &
                        number_4)
                      end if
                      call getjjJMinMax(                               &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il-1,&
                      N3_big_jj,JJTM_min,JJTM_max)
                      call getjjJMinMax(                               &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+1,&
                      N4_big_jj,JJTP_min,JJTP_max)
              do JJTM = JJTM_min, JJTM_max, 2
              do JJTP = JJTP_min, JJTP_max, 2
                      J_2_min =ITREXG(J_1,JJTM,JJTP,iJ_total,IKK)
                      if(IKK.LE.0) cycle
                      J_2_max = J_2_min+IKK-1
              do J_2 = J_2_min, J_2_max, 2
!
                      do iterm_3 = 1,number_3
                        if(jj_3_term(iterm_3)%subshellJ /= JJTM) cycle
                        do iterm_4 = 1,number_4
                          if(jj_4_term(iterm_4)%subshellJ /= JJTP) cycle
                          if(inosubc /= 2) cycle
                          icsf_jj2 = icsf_jj2 + 1
                          expansion_cfg_jj2%csfs(icsf_jj2)%nosubc=     &
                                                              inosubc
                          allocate(expansion_cfg_jj2%csfs(icsf_jj2)%   &
                          subc_cfg(2))
                          allocate(expansion_cfg_jj2%csfs(icsf_jj2)%   &
                          subc(2))
                          allocate(expansion_cfg_jj2%csfs(icsf_jj2)%   &
                          iM1(2))
                          allocate(expansion_cfg_jj2%csfs(icsf_jj2)%   &
                          iM2(2))
                          allocate(expansion_cfg_jj2%csfs(icsf_jj2)%   &
                          iJ(2))
                          expansion_cfg_jj2%csfs(icsf_jj2)%iM1(1)      &
                               = 1000*N1_big_jj + N3_big_jj
                          expansion_cfg_jj2%csfs(icsf_jj2)%iM1(2)      &
                               = 1000*iterm_1 + iterm_3
                          expansion_cfg_jj2%csfs(icsf_jj2)%iM2(1)      &
                               = 1000*N2_big_jj + N4_big_jj
                          expansion_cfg_jj2%csfs(icsf_jj2)%iM2(2)      &
                               = 1000*iterm_2 + iterm_4
                          expansion_cfg_jj2%csfs(icsf_jj2)%iJ(1)=J_1
                          expansion_cfg_jj2%csfs(icsf_jj2)%iJ(2)=J_2
                          expansion_cfg_jj2%csfs(icsf_jj2)%subc_cfg(1)=&
                          expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)
!
                          expansion_cfg_jj2%csfs(icsf_jj2)%subc(1) =   &
                          expansion_cfg_LS%csfs(icsf_nr)%subc(1)
!
                          expansion_cfg_jj2%csfs(icsf_jj2)%subc_cfg(2)=&
                          expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)
!
                          expansion_cfg_jj2%csfs(icsf_jj2)%subc(2) =   &
                          expansion_cfg_LS%csfs(icsf_nr)%subc(2)
                        end do
                      end do
              end do
              end do
              end do
                    end if
                          end do
                        end do
                      end do
                    end if
                 end do
           end do
         end do
!
      end if
!      write(*,*) '      end subroutine form_csfs_jj2'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine  define_number_of_csfs_jj2(irez)
!                                                                      *
!     This subroutine defines the number of csfs in jj2 coupling       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
         implicit none
         integer, intent(out)::irez
!GG!         real(kind=dp) :: coef
         integer :: J_1, J_1_min, J_1_max
         integer :: J_2, J_2_min, J_2_max, J2_summ
         integer :: IKK,JJTM,JJTP,JJTM_min,JJTM_max,JJTP_min,JJTP_max
         integer :: inosubc, isubc, isubc_aviable
         integer :: icsf, icsf_nr, ITTK, ITREXG
         integer :: I, II, N1_MAX, N2_MAX
         integer :: N1_big, N2_big, N3_big, N4_big
         integer :: number_1, number_2, number_3, number_4
         integer :: I_Count1, I_Count2
         integer :: iterm_1, iterm_2, iterm_3, iterm_4
         integer :: N1_big_jj, N2_big_jj, N3_big_jj, N4_big_jj
         integer :: Q_Term1, Nr_Term1, Q_Term2, Nr_Term2
         type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
         type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
!      write(*,*) '      subroutine define_number_of_csfs_jj2'
         isubc_aviable = 2
         if(expansion_cfg_LS%csfs(1)%nosubc .gt. isubc_aviable) then
           write(*,*) 'STOP at subroutine define_number_of_csfs_jj2 ', &
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
              Q_Term1 =                                                &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il+1-&
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
           if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il == 3 .and. &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big>=3) then
            Nr_Term2 = expansion_cfg_LS%csfs(icsf_nr)%subc(2)%inr
            Q_Term2 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il,             &
            expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big,Nr_term2,&
            expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iL,                 &
            expansion_cfg_LS%csfs(icsf_nr)%subc(2)%iS)
           else
              Q_Term2 =                                                &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+1-&
                     expansion_cfg_LS%csfs(icsf_nr)%subc(2)%inr
              Nr_Term2 = 0
           end if
           N2_MAX=2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il
           if (N2_MAX == 0) then
              N3_big = 0
              N4_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big
              I_Count2 = 1
           else if(expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big   &
                                                        <= N2_MAX) then
              N3_big = expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big
              N4_big= 0
              I_Count2 = N3_big + 1
           else
              N3_big=N2_MAX
              N4_big=                                                  &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%iN_big-N3_big
              I_Count2 =                                               &
                2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+2      &
                -N4_big+1
           end if
           call getJMinMax(                                            &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%il,         &
                expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(1)%iN_big,     &
                J_1_min,J_1_max)
           do J_1 = J_1_min, J_1_max, 2
                 do I = 1,I_Count1,1
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
!GG!                     Nr_Term1,                                         &
!GG!                     Q_Term1,                                          &
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
!GG!                          if (coef == ZERO_dp) cycle
!
                 do II = 1,I_Count2,1
                    if (N3_big-II+1 .ge. 0) then
                      if(II == 1) then
                        if (N2_MAX /= 0) then
                           N3_big_jj = N3_big
                           call gettermjj(                      &
                           2*expansion_cfg_LS%csfs(icsf_nr)%    &
                           subc_cfg(2)%il-1,N3_big_jj,jj_3_term,&
                           number_3)
                        else
                           N3_big_jj = 0
                           number_3 = 1
                           jj_3_term(1)%j = 0
                           jj_3_term(1)%Q = 0
                           jj_3_term(1)%subshellJ = 0
                        end if
                        N4_big_jj = N4_big
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il+1,N4_big_jj,jj_4_term,  &
                        number_4)
                      else
                        N3_big_jj = N3_big-II+1
                        N4_big_jj = N4_big+II-1
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il-1,N3_big_jj,jj_3_term,  &
                        number_3)
                        call gettermjj(                        &
                        2*expansion_cfg_LS%csfs(icsf_nr)%      &
                        subc_cfg(2)%il+1,N4_big_jj,jj_4_term,  &
                        number_4)
                      end if
                      call getjjJMinMax(                               &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il-1,&
                      N3_big_jj,JJTM_min,JJTM_max)
                      call getjjJMinMax(                               &
                     2*expansion_cfg_LS%csfs(icsf_nr)%subc_cfg(2)%il+1,&
                      N4_big_jj,JJTP_min,JJTP_max)
              do JJTM = JJTM_min, JJTM_max, 2
              do JJTP = JJTP_min, JJTP_max, 2
                      J_2_min =ITREXG(J_1,JJTM,JJTP,iJ_total,IKK)
                      if(IKK.LE.0) cycle
                      J_2_max = J_2_min+IKK-1
              do J_2 = J_2_min, J_2_max, 2
!
                      do iterm_3 = 1,number_3
                        if(jj_3_term(iterm_3)%subshellJ /= JJTM) cycle
                        do iterm_4 = 1,number_4
                          if(jj_4_term(iterm_4)%subshellJ /= JJTP) cycle
                           irez = irez +1
                        end do
                      end do
              end do
              end do
              end do
                    end if
                          end do
                        end do
                      end do
                    end if
                 end do
           end do
         end do
!      write(*,*) '     end subroutine define_number_of_csfs_jj2'
         end subroutine  define_number_of_csfs_jj2
      end subroutine form_csfs_jj2
!
!***********************************************************************
!                                                                      *
      subroutine matrix_LS_jj2(icsf_LS,icsf_jj2,rez)
!                                                                      *
!     This subroutine calculates the transformation                    *
!     matrix element between two csfs  (in LS and jj2 couplings)       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
      integer,intent(in):: icsf_LS,icsf_jj2
      real(kind=dp),intent(out)::rez
      real(kind=dp) :: rez_1, rez_2, rez_3
      integer :: L1,IS1,L2,IS2,JM2,JP2
      integer :: L,IS,J_1,J_2,J, J_2_min, J_2_max, J_12
      integer :: nosubc, N1, N2, N3, N4
      integer :: number_1, number_2, number_3, number_4
      integer :: num_1, num_2, num_3, num_4
      integer :: Q_Term1, Nr_Term1, Q_Term2, Nr_Term2
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
!      write(*,*) '      subroutine matrix_LS_jj2'
      nosubc=expansion_cfg_LS%csfs(1)%nosubc
      rez = ZERO_dp
      if(equivalent_LSjj2(expansion_cfg_LS%csfs(icsf_LS),            &
         expansion_cfg_jj2%csfs(icsf_jj2))) then
         L1   = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iL
         IS1  = expansion_cfg_LS%csfs(icsf_LS)%subc(1)%iS
         L2   = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iL
         IS2  = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iS
         L    = expansion_cfg_LS%csfs(icsf_LS)%iM1(2)
         iS   = expansion_cfg_LS%csfs(icsf_LS)%iM2(2)
         J_1  = expansion_cfg_jj2%csfs(icsf_jj2)%iJ(1)
         J_12 = expansion_cfg_jj2%csfs(icsf_jj2)%iJ(2)
         J    = states%states(expansion_cfg_LS%nr_of_state)%J

         N1   = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM1(1),2,1000)
         N2   = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM2(1),2,1000)
         number_1 = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM1(2),2,1000)
         number_2 = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM2(2),2,1000)
         N3   = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM1(1),1,1000)
         N4   = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM2(1),1,1000)
         number_3 = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM1(2),1,1000)
         number_4 = JTHN(expansion_cfg_jj2%csfs(icsf_jj2)%iM2(2),1,1000)
         J_2_min = iabs(J_1 - J)
         J_2_max = J_1 + J
         if(J_2_min > J_2_max) return
         if(expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il == 3 .and.   &
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
         if(expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il == 3 .and.   &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%iN_big>=3) then
            Nr_Term2 = expansion_cfg_LS%csfs(icsf_LS)%subc(2)%inr
            Q_Term2 = gettermLSQ(                                      &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il,             &
            expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%iN_big,Nr_term2,&
            expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iL,                 &
            expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iS)
         else
            Q_Term2 =                                                  &
                   2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il+1-  &
                   expansion_cfg_LS%csfs(icsf_LS)%subc(2)%inr
            Nr_Term2 = 0
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
         if (expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il /= 0) then
            call gettermjj                                             &
               (2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il-1,     &
               N3,jj_3_term, num_3)
         else
            jj_3_term(1)%j = 0
            jj_3_term(1)%Q = 0
            jj_3_term(1)%subshellJ = 0
         end if
         call gettermjj                                                &
            (2*expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il+1,        &
            N4,jj_4_term,num_4)
!
         rez_1 = coefLSjj(                                         &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%il,    &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(1)%iN_big,&
                 Nr_Term1,                                         &
                 Q_Term1,                                          &
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
         if (rez_1 == ZERO_dp) return
         JM2 = jj_3_term(number_3)%subshellJ
         JP2 = jj_4_term(number_4)%subshellJ
         do J_2 = J_2_min, J_2_max, 2
            rez_2 = coefLSjj(                                      &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%il,    &
                 expansion_cfg_LS%csfs(icsf_LS)%subc_cfg(2)%iN_big,&
                 Nr_Term2,                                         &
                 Q_Term2,                                          &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iL,        &
                 expansion_cfg_LS%csfs(icsf_LS)%subc(2)%iS,        &
                 J_2,                                              &
                 jj_3_term(number_3)%j,                            &
                 N3,                                               &
                 jj_3_term(number_3)%Q,                            &
                 jj_3_term(number_3)%subshellJ,                    &
                 jj_4_term(number_4)%j,                            &
                 jj_4_term(number_4)%Q,                            &
                 jj_4_term(number_4)%subshellJ)
            if (rez_2 == ZERO_dp) cycle
            call jj2PER(L1,IS1,L2,IS2,L,IS,J_1,J_2,J,JM2,JP2,J_12,rez_3)
            rez = rez + rez_1 * rez_2 * rez_3
         end do
      end if
!      write(*,*) '      end subroutine matrix_LS_jj2'
      end subroutine matrix_LS_jj2
!
!***********************************************************************
!                                                                      *
      subroutine dotransform_LSjj2
!                                                                      *
!     This subroutine calculates the weights of the                    *
!     expansions in the jj2 coupling                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp)::coeff_jj2,coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_jj2
!      write(*,*) '      subroutine dotransform_LSjj2'
      melement=0.
      sum_of_state = 0.
      do icsf_jj2=1,expansion_cfg_jj2%size
         coeff_jj2 = 0.
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_jj2(icsf_LS,icsf_jj2,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_jj2 = coeff_jj2 + coeff_LS*melement
         end do
         if((dabs(coeff_jj2)-dp_coeff_precision).gt.ONE_dp) then
            write(*,*)'possible error at subroutine dotransform_LSjj2:',&
               ' coeff_jj2=',coeff_jj2
            write(iwrite_log,*)'possible error at subroutine ',        &
               'dotransform_LSjj2: coeff_jj2=',coeff_jj2
         end if
         expansion_cfg_jj2%coeffs(icsf_jj2)= coeff_jj2
         sum_of_state = sum_of_state + coeff_jj2*coeff_jj2
      end do
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_LSjj2: ',  &
            'sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine ',           &
            'dotransform_LSjj2: sum_of_state=',sum_of_state
      end if
!      write(*,*) '      end subroutine dotransform_LSjj2'
      end subroutine dotransform_LSjj2
!
!***********************************************************************
!                                                                      *
      subroutine print_cfg_LSjj2(itype)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the unit "iwrite_cfg_expansions_jj2"                *
!     (see module "Coupling_constants")                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in) :: itype
      integer             :: icsf, isubc, j!, JTHN
      integer             :: num_1, num_2, num_3, num_4
      integer             :: number_1, number_2, number_3, number_4
      character(len=1)    :: CVAL
      character(len=2)    :: LSJVAL
      character(len=4)    :: JVAL
      type(subshell_term), dimension(1:63) :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63) :: jj_3_term, jj_4_term
!
!     LS - Coupling
!
!      write(*,*) '      subroutine print_cfg_LSjj2'
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
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_jj123,                &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
                 'To many couled shells in Coupling_transform_cfg_LSjj2'
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
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_jj123,               &
                     '(12x,7(1x,i2,a1,i1))')                         &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,      &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL),&
                      expansion_cfg_LS%csfs(icsf)%subc(j)%inr,       &
                      j=1,2),                                        &
                      expansion_cfg_LS%csfs(icsf)%iM2(2)+1,          &
                      CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(2))
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
!     jj2 - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_jj123,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_jj123,*) &
         'state Nr.',expansion_cfg_jj2%nr_of_state
         write(iwrite_cfg_expansions_jj123,'(3x,a3,a4,a9,3x,f15.8)')   &
         'J =', JVAL(states%states(expansion_cfg_jj2%nr_of_state)%J),  &
         ' Energy =',states%states(expansion_cfg_jj2%nr_of_state)%energy
         write(iwrite_cfg_expansions_jj123,'(3x,a29,i2)')              &
         'expansion size (jj2 coupling): ', expansion_cfg_jj2%size
         write(iwrite_cfg_expansions_jj123,*)''
         if(associated(expansion_cfg_jj2%csfs)) then
            do icsf=1,expansion_cfg_jj2%size
               if(associated(expansion_cfg_jj2%csfs(icsf)%subc_cfg)) then
                  if (expansion_cfg_jj2%csfs(icsf)%nosubc == 2) then
                     number_1 =                                        &
                       JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(2),2,1000)
                     number_2 =                                        &
                       JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(2),2,1000)
                     number_3 =                                        &
                       JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(2),1,1000)
                     number_4 =                                        &
                       JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(2),1,1000)
                     if(expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%il/=0)&
                                                                   then
                        call gettermjj                                 &
                        (2*expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%il &
                        -1,                                            &
                      JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(1),2,1000),&
                        jj_1_term,num_1)
                     else
                        jj_1_term(number_1)%j = -1
                        jj_1_term(number_1)%subshellJ = 0
                        jj_1_term(number_1)%nu  = 0
                     end if
                     call gettermjj                                    &
                     (2*expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%il+1, &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(1),2,1000), &
                     jj_2_term,num_2)
!
                     if(expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%il/=0)&
                                                                   then
                        call gettermjj                                 &
                        (2*expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%il &
                        -1,                                            &
                      JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(1),1,1000),&
                        jj_3_term,num_3)
                     else
                        jj_3_term(number_3)%j = -1
                        jj_3_term(number_3)%subshellJ = 0
                        jj_3_term(number_3)%nu  = 0
                     end if
                     call gettermjj                                    &
                     (2*expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%il+1, &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(1),1,1000), &
                     jj_4_term,num_4)
                     write(iwrite_cfg_expansions_jj123,                &
                     '(1x,i5,3x,i2,a2,"(",i1,")",                      &
                     3(i2,a2,"(",i1,")",6x),a6,f10.7)')                &
                     icsf,                                             &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%il,      &
                     jj_1_term(number_1)%j),                           &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(1),2,1000), &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(1)%il,      &
                     jj_2_term(number_2)%j),                           &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(1),2,1000), &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%il,      &
                     jj_3_term(number_3)%j),                           &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM1(1),1,1000), &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%in,      &
                     LSJVAL(                                           &
                     expansion_cfg_jj2%csfs(icsf)%subc_cfg(2)%il,      &
                     jj_4_term(number_4)%j),                           &
                     JTHN(expansion_cfg_jj2%csfs(icsf)%iM2(1),1,1000), &
                     'coeff:',expansion_cfg_jj2%coeffs(icsf)
                  end if
               else
                  write(iwrite_cfg_expansions_jj123,'(9x,a51)') &
                  'expansion_cfg_jj2%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_jj2%csfs(icsf)%subc)     &
                  .and.associated(expansion_cfg_jj2%csfs(icsf)%iM1) &
                  .and.associated(expansion_cfg_jj2%csfs(icsf)%iM2)) then
                  if (expansion_cfg_jj2%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_jj123,                &
                     '(9x,2(a4,1x,i1,1x),"[",a4,"]",                   &
                          2(a4,1x,i1,1x,"[",a4,"]"))')                 &
                     JVAL(jj_1_term(number_1)%subshellJ),              &
                     jj_1_term(number_1)%nu,                           &
                     JVAL(jj_2_term(number_2)%subshellJ),              &
                     jj_2_term(number_2)%nu,                           &
                     JVAL(expansion_cfg_jj2%csfs(icsf)%iJ(1)),         &
                     JVAL(jj_3_term(number_3)%subshellJ),              &
                     jj_3_term(number_3)%nu,                           &
                     JVAL(expansion_cfg_jj2%csfs(icsf)%iJ(2)),         &
                     JVAL(jj_4_term(number_4)%subshellJ),              &
                     jj_4_term(number_4)%nu,                           &
                    JVAL(states%states(expansion_cfg_jj2%nr_of_state)%J)
                  end if
               else
                  write(iwrite_cfg_expansions_jj123,'(9x,a63)') &
       'expansion_cfg_jj2%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_jj123,*)' '
            end do
         else
            write(iwrite_cfg_expansions_jj123,*) &
            'expansion_cfg_jj2%csfs NOT associated'
         end if
            write(iwrite_cfg_expansions_jj123,*) &
            '-------------------------------------'
      end if
!      write(*,*) '      end subroutine print_cfg_LSjj2'
      end subroutine print_cfg_LSjj2
!
      end module Coupling_transform_cfg_LSjj2
