!
!***********************************************************************
!                                                                      *
      module Coupling_evaluation
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
      use Coupling_inside_shell
!
      public  :: main_evaluation
      private :: generate_classification_data
      private :: print_classification_data
      private :: define_nr_of_max_csf
      private :: generate_evaluation_data
      private :: print_evaluation_data
      private :: generate_J_structure
      private :: delete_J_structure
      private :: generate_suggestions
      public  :: getchLS
      public  :: getchjj
      public  :: getchXLS
      private :: CONVRT_DOUBLE
      private :: CONVRT
      public  :: Seniority
      public  :: Seniority_JJ
!
!--- definition of types   ------------------------
!
      type::J_list
         integer::J !number of states with given J
         type(list)::states !i=1..
      end type J_list
!
      type:: J_lists
         integer::nr_of_Js
         type(J_list),dimension(:),pointer::Js_states !i=1..nr_of_Js
      end type J_lists
!
!--- definition of variables   ------------------------
!
      type(J_lists), private :: J_structure
      character(len=1),dimension(0:20),private, parameter :: L_string =&
               (/ "S","P","D","F","G","H","I","K","L","M","N","O","Q", &
                  "R","T","U","V","W","X","Y","Z" /)
      character(len=1),dimension(0:20),private, parameter :: ll_string=&
               (/ "s","p","d","f","g","h","i","k","l","m","n","o","q", &
                  "r","t","u","v","w","x","y","z" /)
! (it contains the lists of serial numbers of states
! with certain J values)

contains

      subroutine main_evaluation(nr_of_j_count,nr_of_j,                &
                      iwrite_cl,iwrite_cl_eval,iwrite_sugg,print_level)
!--------------------------------------------------------------------
! This is managerial subroutine
!--------------------------------------------------------------------
      implicit none
!integer::error
!integer::itype
      integer, intent(in) :: nr_of_j_count,nr_of_j
      integer, intent(in) :: iwrite_cl,iwrite_cl_eval
      integer, intent(in) :: iwrite_sugg, print_level
!
!      write(*,*) '  subroutine main_evaluation'
      write(iwrite_log,*)                                              &
      '------------------------------------------------------'
      call generate_classification_data
      if(print_level.gt.0) call print_classification_data(iwrite_cl)
      call generate_J_structure
      call generate_evaluation_data(nr_of_j_count,nr_of_j)
      if(print_level.gt.0) call print_evaluation_data(iwrite_cl_eval)
      write(iwrite_log,*)                                              &
      '------------------------------------------------------'
      if(print_level > 2) call generate_suggestions(iwrite_sugg)
      call delete_J_structure
!      write(*,*) '  end subroutine main_evaluation'
!      write(iwrite_log,*)""
!
      end subroutine main_evaluation
!
!-------------------------------------------------------------------
!
      subroutine generate_suggestions(iwrite)
!--------------------------------------------------------------------
! This subroutine generates the suggestion of the optimal spectrum
! classification and prints it to the file
!--------------------------------------------------------------------
      implicit none
      integer, intent(in):: iwrite
!integer::error
      integer::icoupling, recomended_coupling
!
!      write(*,*) '    subroutine generate_suggestions'
      if(all_evaluations%nr_of_couplings .ne.                          &
                             all_classifications%nr_of_couplings) then
         write(*,*) 'STOP at subroutine generate_suggestions:', &
         ' all_evaluations%nr_of_couplings .ne.',               &
         ' all_classifications%nr_of_couplings'
         stop
      end if
!
!find the best coupling
!
      if(is_there_one_to_one_couplings()) then
         write(iwrite_log,*)'There is one-to-one coupling'
         recomended_coupling=first_one_to_one()
         do icoupling=1, all_evaluations%nr_of_couplings,1
            if(all_classifications%couplings(icoupling)%is_one_to_one) then
               if(all_evaluations%R(icoupling)%R_value .lt.           &
                  all_evaluations%R(recomended_coupling)%R_value)     &
                                         recomended_coupling=icoupling
            end if
         end do
      else
         write(iwrite_log,*)'There is NOT one-to-one coupling !!!'
         write(*,*)' There is NOT one-to-one coupling'
         recomended_coupling = 1
         do icoupling=1, all_evaluations%nr_of_couplings,1
            if(all_evaluations%R(icoupling)%R_value .lt.              &
               all_evaluations%R(recomended_coupling)%R_value)        &
                                         recomended_coupling=icoupling
         end do
      end if
      call print_suggestions(iwrite)
!      write(*,*) '    end subroutine generate_suggestions'
!
      contains
!
         function is_there_one_to_one_couplings()      result(rez)
!--------------------------------------------------------------------
! This function defines whether is there one-to-one
! couplings or not
!--------------------------------------------------------------------
         integer :: icoupling
         logical :: rez
         rez=.false.
         do icoupling=1, all_evaluations%nr_of_couplings, 1
            if(all_classifications%couplings(icoupling)%is_one_to_one) &
                                                                    then
               rez=.true.
               exit
            end if
         end do
         end function is_there_one_to_one_couplings
!
         function first_one_to_one()      result(irez)
!--------------------------------------------------------------------
! This function defines the smallest serial number
! of one-to-one	coupling
!--------------------------------------------------------------------
         integer :: icoupling
         integer :: irez
         do icoupling=1, all_evaluations%nr_of_couplings, 1
            if(all_classifications%couplings(icoupling)%is_one_to_one) &
                                                                    then
               irez=icoupling
               exit
            end if
         end do
         end function first_one_to_one
!
         subroutine print_suggestions(iwrite)
!--------------------------------------------------------------------
! This subroutine prints suggestions to the unit "iwrite"
!--------------------------------------------------------------------
         integer, intent(in):: iwrite
!
         write(*,*)''
         write(*,*)'We suggest to classify spectra in: ',              &
            coupling_descriptions(all_expansions%coupling_expansions(  &
            recomended_coupling)%icoupling)%long_name
         write(iwrite,*)                                               &
         '************************************************'
         write(iwrite,*)'We suggest to classify spectra in: ',         &
            coupling_descriptions(all_expansions%coupling_expansions(  &
            recomended_coupling)%icoupling)%long_name
         write(iwrite,*)' '
         write(iwrite,*)'('
!
         if(all_classifications%couplings(recomended_coupling)%        &
            is_one_to_one) then
            write(iwrite,*) '   coupling is one-to-one'
         else
            write(iwrite,*) '   coupling is NOT one-to-one'
         end if
         write(iwrite,'(3x,a4,f9.6)')'R = ',                           &
                      all_evaluations%R(recomended_coupling)%R_value
         write(iwrite,*)')'
         write(iwrite,*)'Corresponding classification presented below.'
         write(iwrite,*)                                               &
         '************************************************'
!
         call print_classification(iwrite, recomended_coupling)
!
         end subroutine print_suggestions
      end subroutine generate_suggestions
!
!-------------------------------------------------------------------
!
      subroutine generate_evaluation_data(nr_of_j_count,nr_of_j)
!--------------------------------------------------------------------
! This subroutine evaluates all the classifications. The subroutine
! calculates all the R_values and P_values for all the couplings.
!--------------------------------------------------------------------
!
      implicit none
      integer, intent(in) :: nr_of_j_count,nr_of_j
      integer::icoupling, istate, iJ, i2Jp1, isum2Jp1p
      integer, save :: Total_States
      real(kind=dp)::RJ_tmp
!
!      write(*,*) '    subroutine generate_evaluation_data'
!
      all_evaluations%nr_of_couplings =                                &
                     all_expansions%nr_of_coupling_expansions
      if(nr_of_j_count == 1) then
         allocate(all_evaluations%R(all_evaluations%nr_of_couplings))
         Total_States = 0
      end if
      do icoupling=1, all_evaluations%nr_of_couplings
         all_evaluations%R(icoupling)%nr_of_J=J_structure%nr_of_Js
!GG         if(nr_of_j_count == 1) then
!GG            allocate(all_evaluations%R(icoupling)%RJ_values(           &
!GG                                  all_evaluations%R(icoupling)%nr_of_J))
!GG            allocate(all_evaluations%R(icoupling)%J(all_evaluations%   &
!GG                                                  R(icoupling)%nr_of_J))
!GG         end if
!calculate RJ values, ...
!
         isum2Jp1p=0
!GG         all_evaluations%R(icoupling)%R_value=ZERO_dp
         if(nr_of_j_count == 1) all_evaluations%R(icoupling)%P_value=ZERO_dp
!
!GG         do iJ=1, all_evaluations%R(icoupling)%nr_of_J
!GG            i2Jp1=J_structure%Js_states(iJ)%J  + 1
!GG            all_evaluations%R(icoupling)%J(iJ) =                       &
!GG                                         J_structure%Js_states(iJ)%J
!GG            if(J_structure%Js_states(iJ)%states%list_size.gt.1) then
!GG               isum2Jp1p = isum2Jp1p + i2Jp1
!GG               all_evaluations%R(icoupling)%RJ_values(iJ) = ZERO_dp
!
!GG               do istate=1, J_structure%Js_states(iJ)%states%list_size
!Here calculation according to Rudzikas-Ciplys (3.8)
!will be performed
! here dp_coeff_precision_t10 added to escape the case of
! negative argument of dsqrt()
!
!GG                  RJ_tmp=dsqrt(TWO_dp*(ONE_dp+dp_coeff_precision_t10 - &
!GG                     dabs(all_classifications%couplings(icoupling)%    &
!GG                     states(J_structure%Js_states(iJ)%states%items(    &
!GG                     istate))%coeff)))
!GG                  all_evaluations%R(icoupling)%RJ_values(iJ) =         &
!GG                     all_evaluations%R(icoupling)%RJ_values(iJ) + RJ_tmp
!
!GG               end do !istate
!GG               all_evaluations%R(icoupling)%RJ_values(iJ) =            &
!GG                  all_evaluations%R(icoupling)%RJ_values(iJ) /         &
!GG                  DBLE(J_structure%Js_states(iJ)%states%list_size)
!
!GG               all_evaluations%R(icoupling)%R_value =                  &
!GG                  all_evaluations%R(icoupling)%R_value +               &
!GG                  all_evaluations%R(icoupling)%RJ_values(iJ)*DBLE(i2Jp1)
!
!GG            else ! ...states%list_size.gt.1
!
!GG               if(J_structure%Js_states(iJ)%states%list_size.gt.0) then
!GG                  isum2Jp1p = isum2Jp1p + i2Jp1
!GG                  all_evaluations%R(icoupling)%RJ_values(iJ)= ZERO_dp
!all_evaluations%R(icoupling)%RJ_values(iJ)= ONE_dp
!GG               else
!GG                  all_evaluations%R(icoupling)%RJ_values(iJ)=ZERO_dp
!GG                  all_evaluations%R(icoupling)%J(iJ)=0
!GG                  write(*,*) 'possible error in subroutine ',          &
!GG                     'generate_evaluation_data: ',                     &
!GG                  '.not.J_structure%Js_states(iJ)%states%list_size.gt.0'
!GG                  write(iwrite_log,*) 'possible error in ',            &
!GG                     'subroutine generate_evaluation_data: ',          &
!GG                  '.not.J_structure%Js_states(iJ)%states%list_size.gt.0'
!GG               end if !  ...states%list_size.gt.0
!
!GG            end if ! ...states%list_size.gt.1
!
!GG         end do ! iJ
!
         all_evaluations%R(icoupling)%R_value =                        &
                 all_evaluations%R(icoupling)%R_value / DBLE(isum2Jp1p)
!
!find out P_value
!ne pats optimaliausias budas bet kolkas tiks ...
         if(all_classifications%couplings(icoupling)%nr_of_states >    &
                                                                 0) then
            if(icoupling == 1) Total_States = Total_States +           &
                          all_classifications%couplings(1)%nr_of_states
            do istate=1,all_classifications%couplings(icoupling)%      &
                                                           nr_of_states
               all_evaluations%R(icoupling)%P_value = all_evaluations% &
                  R(icoupling)%P_value + all_classifications%couplings(&
                  icoupling)%states(istate)%coeff*all_classifications% &
                  couplings(icoupling)%states(istate)%coeff
            end do !istate
            if(nr_of_j_count == nr_of_j) then
               all_evaluations%R(icoupling)%P_value = all_evaluations% &
                  R(icoupling)%P_value/DBLE(Total_States)
            end if
         else
            all_evaluations%R(icoupling)%P_value=ONE_dp
         end if
!
      end do ! icoupling
!
!      write(*,*) '    end subroutine generate_evaluation_data'
!
      end subroutine generate_evaluation_data
!
!-------------------------------------------------------------------
!
      subroutine print_evaluation_data(iwrite)
!--------------------------------------------------------------------
! This subroutine prints the evaluation data to the unit "iwrite"
!--------------------------------------------------------------------
!
      implicit none
      integer, intent(in):: iwrite
      integer::icoupling, iJ
      character(len=4) :: JVAL
!real(kind=dp)::coeff
!integer::itype
!
!      write(*,*) '    subroutine print_evaluation_data'
      write(iwrite,*) 'all_evaluations%nr_of_couplings', &
                       all_evaluations%nr_of_couplings
      if(associated(all_evaluations%R)) then
         do icoupling=1, all_evaluations%nr_of_couplings
            write(iwrite, *)'----------------------------------------'
            write(iwrite, *)                                           &
            coupling_descriptions(all_expansions%                      &
                    coupling_expansions(icoupling)%icoupling)%long_name
            write(iwrite, *)'  nr_of_J:',                              &
                                  all_evaluations%R(icoupling)%nr_of_J
            write(iwrite,*)' '
!GG            do iJ=1, all_evaluations%R(icoupling)%nr_of_J
!GG               write(iwrite,'(4x,a7,a4,f9.5)')'J, RJ: ',               &
!GG                             JVAL(all_evaluations%R(icoupling)%J(iJ)), &
!GG                           all_evaluations%R(icoupling)%RJ_values(iJ)
!GG               write(iwrite,*)' '
!GG            end do
!GG            write(iwrite,'(2x,a3,f9.5)')'R: ',                         &
!GG                                 all_evaluations%R(icoupling)%R_value
            write(iwrite,*)' '
         end do
      end if
!      write(*,*) '    end subroutine print_evaluation_data'
!
      end subroutine print_evaluation_data
!
!-------------------------------------------------------------------
!
      subroutine generate_J_structure
!--------------------------------------------------------------------
! This subroutine generates the structure "J_structure"
! (it contains the lists of serial numbers of states
! with certain J values)
!--------------------------------------------------------------------
      implicit none
!integer::error
      integer:: istate !, icoupling ,istate2 !, imax_csf_nr
      integer::iJ, nr_of_iJ, inr_of_states_of_iJ
!real(kind=dp)::coeff
!integer::itype
!
!      write(*,*) '    subroutine generate_J_structure'
!
      J_structure%nr_of_Js=count_nr_of_Js()
      allocate(J_structure%Js_states(J_structure%nr_of_Js))
      if(associated(states%states).and.states%nr_of_states.gt.0) then
         nr_of_iJ=1
         J_structure%Js_states(nr_of_iJ)%J=states%states(1)%J
         inr_of_states_of_iJ=1
         do istate=2, states%nr_of_states,1
            if(states%states(istate)%J .eq.                            &
                                 J_structure%Js_states(nr_of_iJ)%J) then
               inr_of_states_of_iJ=inr_of_states_of_iJ+1
            else
!fix the results of count
               J_structure%Js_states(nr_of_iJ)%states%list_size =      &
                                                     inr_of_states_of_iJ
               allocate(J_structure%Js_states(nr_of_iJ)%states%        &
                items(J_structure%Js_states(nr_of_iJ)%states%list_size))
!start new count for nr_of_iJ=nr_of_iJ+1
               nr_of_iJ=nr_of_iJ+1
               if(nr_of_iJ.gt.J_structure%nr_of_Js) then
                  write(*,*)'STOP at subroutine generate_J_structure:',&
                            ' nr_of_iJ.gt.J_structure%nr_of_Js'
                  stop
               end if
               J_structure%Js_states(nr_of_iJ)%J=states%states(istate)%J
               inr_of_states_of_iJ=1
            end if
         end do !istates
!fix the results of last count
         J_structure%Js_states(nr_of_iJ)%states%list_size =            &
                                                     inr_of_states_of_iJ
         allocate(J_structure%Js_states(nr_of_iJ)%states%              &
                items(J_structure%Js_states(nr_of_iJ)%states%list_size))
!fill up this structure
!
         nr_of_iJ = 0
         do iJ=1, J_structure%nr_of_Js
            do istate=1, J_structure%Js_states(iJ)%states%list_size
               nr_of_iJ = nr_of_iJ+1
               J_structure%Js_states(iJ)%states%items(istate)=nr_of_iJ
            end do
         end do
      end if !associated(states%states).and.states%nr_of_states.gt.0
!call print_J_structure
!
!      write(*,*) '    end subroutine generate_J_structure'
!
      contains
!
         function count_nr_of_Js()    result(irez)
!--------------------------------------------------------------------
! This function	counts the number of different J values
!--------------------------------------------------------------------
         implicit none
         integer::irez
         integer::istate, iJtemp
!
         if(associated(states%states).and.states%nr_of_states.gt.0) then
            iJtemp=states%states(1)%J
            irez=1
            do istate=2, states%nr_of_states
               if(states%states(istate)%J.ne.iJtemp) then
                  irez=irez+1
                  iJtemp=states%states(istate)%J
               end if
            end do !istates
         end if
         end function count_nr_of_Js
!
!subroutine print_J_structure
!	integer::istate, iJ
!	write(*,*)'   -------  print J_structure  ----------'
!	write(*,*) '     J_structure%nr_of_Js : ', J_structure%nr_of_Js
!	do iJ=1, J_structure%nr_of_Js
!		write(*,*) '     J: ', J_structure%Js_states(iJ)%J
!		write(*,*) '       nr_of_states: ', J_structure%Js_states(iJ)%states%list_size
!		do istate=1, J_structure%Js_states(iJ)%states%list_size
!			write(*,*) '          istate: ', J_structure%Js_states(iJ)%states%items(istate)
!		end do
!	end do
!	write(*,*)'   -------  end print J_structure  ----------'
!end subroutine print_J_structure
!
      end subroutine generate_J_structure
!
!-------------------------------------------------------------------
!
      subroutine delete_J_structure
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of "J_structure"
!--------------------------------------------------------------------
      implicit none
      integer::error
      integer:: iJ
!      write(*,*) '    subroutine delete_J_structure'
      if(associated(J_structure%Js_states)) then
         do iJ=1, J_structure%nr_of_Js
            if(associated(J_structure%Js_states(iJ)%states%items))     &
           deallocate(J_structure%Js_states(iJ)%states%items,STAT=error)
         end do
         deallocate(J_structure%Js_states, STAT=error)
      end if
!      write(*,*) '    end subroutine delete_J_structure'
!
      end subroutine delete_J_structure
!
!-----------------------------------------------------------------
!
      subroutine generate_classification_data
!--------------------------------------------------------------------
! This subroutine generates the classification data
! "all_classifications".
! In this version of the program the states is
! classified by csf with the biggest weight
! in the expansion
!--------------------------------------------------------------------
      implicit none
      integer :: icoupling, istate, istate2, imax_csf_nr, indent
      real(kind=dp):: coeff
!
!      write(*,*) '    subroutine generate_classification_data'
!
      all_classifications%nr_of_couplings =                            &
                               all_expansions%nr_of_coupling_expansions
      allocate(all_classifications%couplings(all_classifications%      &
                                                       nr_of_couplings))
!
      do icoupling=1, all_classifications%nr_of_couplings, 1
         indent = 0
         all_classifications%couplings(icoupling)%nr_of_states =       &
         all_expansions%coupling_expansions(icoupling)%nr_of_expansions
         allocate(all_classifications%couplings(icoupling)%            &
          states(all_classifications%couplings(icoupling)%nr_of_states))
!
         do istate=1,                                                  &
                all_classifications%couplings(icoupling)%nr_of_states, 1
          call define_nr_of_max_csf(icoupling,istate,imax_csf_nr, coeff)
!
!unsing this I am SURE that in "b=a" "b" is "empty" csf ...

          if(imax_csf_nr>all_expansions%coupling_expansions(icoupling)&
                                  %expansions(istate)%size) THEN
             print*, "ERROR in generate_classification_datai A"
             print*, "icoupling",icoupling
             print*,"imax_csf_nr=",imax_csf_nr,"istate=",istate
             print*, "all_expansions%coupling_expansions(&
                      icoupling)%expansions(istate)%size", &
            all_expansions%coupling_expansions(icoupling)%expansions(istate)%size
             STOP
          END IF

          all_classifications%couplings(icoupling)%states(istate)%csf =&
          all_expansions%coupling_expansions(icoupling)%               &
                                    expansions(istate)%csfs(imax_csf_nr)
!write(*,'(a32,i2,2x,i2,2x,i2)') 'icoupling, istate, imax_csf_nr: ', icoupling, istate, imax_csf_nr
          all_classifications%couplings(icoupling)%states(istate)%coeff&
          = coeff
!all_classifications%couplings(icoupling)%states(istate)%coeff=all_expansions%coupling_expansions(icoupling)%expansions(istate)%coeffs(imax_csf_nr)
!
         end do
!
         all_classifications%couplings(icoupling)%is_one_to_one = .true.
         if(all_classifications%couplings(icoupling)%nr_of_states .gt. &
                                                                1) then
            do istate=1,(all_classifications%couplings(icoupling)%     &
                                                       nr_of_states-1),1
               do istate2=(istate+1),                                  &
                 all_classifications%couplings(icoupling)%nr_of_states,1
!
                  if(all_classifications%couplings(icoupling)%         &
                                    states(istate)%csf%nosubc == 1) then
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(1)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(1)) cycle
                  else if(all_classifications%couplings(icoupling)%    &
                                    states(istate)%csf%nosubc == 2) then
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(1)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(1)) cycle
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(2)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(2)) cycle
                  else if(all_classifications%couplings(icoupling)%    &
                                    states(istate)%csf%nosubc == 3) then
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(1)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(1)) cycle
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(2)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(2)) cycle
                    if(all_classifications%couplings(icoupling)%states(&
                      istate)%csf%iJ(3)/=all_classifications%couplings(&
                             icoupling)%states(istate2)%csf%iJ(3)) cycle
                  end if
                  if(all_classifications%couplings(icoupling)%states(  &
                     istate)%csf==all_classifications%couplings(       &
                                    icoupling)%states(istate2)%csf) then
                     if(states%states(istate)%J .eq. states%states(    &
                                                        istate2)%J) then
                        all_classifications%couplings(icoupling)%      &
                                                 is_one_to_one = .false.
!GG Issigimimas
                        if(indent == 0) then
                           indent = 1
                           write(iwrite_log,*)""
                           write(iwrite_log,*)            &
                           "There are the states with same ",          &
                           "labeling in ",                      &
                           coupling_descriptions(all_expansions%       &
                           coupling_expansions(icoupling)%icoupling)%  &
                           long_name
                           write(iwrite_log,*)                         &
                           "   States", istate, " and", istate2
                        else
                           write(iwrite_log,*)                         &
                           "         ", istate, "    ", istate2
                        end if
                        exit
                     end if
                  end if
!
               end do
            end do
         end if
         if(indent == 0) then
            write(*,*)"There is one-to-one classification for ",       &
            coupling_descriptions(all_expansions%coupling_expansions(  &
            icoupling)%icoupling)%long_name
         end if
!        print*, "generate_classification_data icoupling pabaiga"
      end do
!
      write(*,*) '    end subroutine generate_classification_data'
!
      end subroutine generate_classification_data
!
!-----------------------------------------------------------------------------
!
      subroutine define_nr_of_max_csf(icoupling,istate,imax_csf_nr,    &
                                                              max_coeff)
!--------------------------------------------------------------------
! This subroutine defines the serial number	of
! the csf (and the corresponding weight coefficient)
! having the biggest weight coefficient in
! the expansion of state "istate" in the coupling "icoupling"
!--------------------------------------------------------------------
      implicit none
      integer, intent(in) :: icoupling, istate
      integer, intent(out) :: imax_csf_nr
      real(kind=dp), intent(out):: max_coeff
      integer::icsf
      imax_csf_nr = 1
      max_coeff = dabs(all_expansions%coupling_expansions(icoupling)%  &
         expansions(istate)%coeffs(1))
      do icsf=2,all_expansions%coupling_expansions(icoupling)%         &
         expansions(istate)%size,1
         if(dabs(all_expansions%coupling_expansions(icoupling)%        &
                    expansions(istate)%coeffs(icsf)).gt.max_coeff) then
            imax_csf_nr = icsf
            max_coeff = dabs(all_expansions%coupling_expansions(       &
               icoupling)%expansions(istate)%coeffs(icsf))
         end if
      end do
      max_coeff= all_expansions%coupling_expansions(icoupling)%        &
         expansions(istate)%coeffs(imax_csf_nr)
      end subroutine define_nr_of_max_csf
!
!------------------------------------------------------------------------
!
      subroutine print_classification_data(iwrite)
!--------------------------------------------------------------------
! This subroutine prints classification data
! to the unit "iwrite"
!--------------------------------------------------------------------
      implicit none
      integer, intent(in):: iwrite
      integer::icoupling, istate, isubc, j
      integer:: iN_big, in, il, Li, Si, inr, L_i, S_i
      integer          :: num_1, num_2, number_1, number_2
      integer          :: num_3, num_4, number_3, number_4
      integer          :: num_5, num_6, number_5, number_6
      character(len=1) :: CVAL
      character(len=2) :: LSJVAL
      character(len=4) :: JVAL
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
      type(subshell_term), dimension(1:63)  :: jj_5_term, jj_6_term
!
!      write(*,*) '    subroutine print_classification_data'
!
      if(associated(all_classifications%couplings)) then
         do icoupling=1, all_classifications%nr_of_couplings
            write(iwrite,*)' '
            write(iwrite,*)     &
            ' ---------------------------------------------- '
            write(iwrite,*)' '
            write(iwrite,*) '  ',coupling_descriptions(all_expansions% &
                    coupling_expansions(icoupling)%icoupling)%long_name
            write(iwrite,*)' '
            if(all_classifications%couplings(icoupling)%is_one_to_one) then
               write(iwrite,*)'   classification is one-to-one'
            else
               write(iwrite,*)'   classification is NOT one-to-one'
            end if
            write(iwrite,*)' '
            if(associated(all_classifications%couplings(icoupling)%    &
                                                          states)) then
               do istate=1,all_classifications%couplings(icoupling)%   &
                                                          nr_of_states,1
                  write(iwrite,*)'   ----------------------------'
                  write(iwrite,*)' '
                  write(iwrite,*)'   state: ', istate
                  write(iwrite,'(5x,a3,a4,a9,3x,f15.8)')   &
                  'J = ', JVAL(states%states(istate)%J),   &
                  ' Energy =',states%states(istate)%energy
                  write(iwrite,*)' '
                  write(iwrite,'(8x,a13)')'configuration'
                  if(associated(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%subc_cfg)) then
                     if(coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj1 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(1),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(1),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM1(1),1,1000),jj_2_term,num_2)
                           write(iwrite,'(9x,i2,a2,"(",i1,")",         &
                           i2,a2,"(",i1,")",6x,a6,f10.7)')             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else if(all_classifications%couplings(icoupling)% &
                                    states(istate)%csf%nosubc == 1) then
                        write(iwrite,                                  &
                           '(6x,3x,i2,a1,"(",i1,")",6x,a6,f10.7)')     &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,CVAL(     &
                           1,all_classifications%couplings(icoupling)% &
                           states(istate)%csf%subc_cfg(1)%il),         &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%iN_big,      &
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LScjj coupling') then
                           number_1 = all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM1(2)
                           number_2 = all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),jj_2_term,num_2)
                           write(iwrite,                               &
                           '(9x,i2,a1,"(",i1,")",2x,             &
                           2(i2,a2,"(",i1,")"),6x,a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           CVAL(1,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il),all_classifications%couplings(          &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           iN_big,                                     &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(2)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_1_term(number_1)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(2)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_2_term(number_2)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),                  &
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj2 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),2,1000)
                           number_3 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),1,1000)
                           number_4 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),2,1000),jj_2_term,num_2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),1,1000),jj_3_term,num_3)
                           else
                              jj_3_term(number_1)%j = -1
                              jj_3_term(number_1)%subshellJ = 0
                              jj_3_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),1,1000),jj_4_term,num_4)
                           write(iwrite,'(9x,i2,a2,"(",i1,")",         &
                           3(i2,a2,"(",i1,")",6x),a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_3_term(number_3)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_4_term(number_4)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff

                     else if(all_classifications%couplings(icoupling)% &
                                    states(istate)%csf%nosubc == 2) then
                        write(iwrite,                                  &
                           '(6x,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')  &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc_cfg(j)%in,CVAL(1,   &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(j)%il),         &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(j)%iN_big,j=1,  &
                           2),'coeff:',all_classifications%couplings(  &
                           icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LS3 coupling' .or.    &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LS coupling' .or.     &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LSJ3 coupling' .or.   &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LK3 coupling' .or.    &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'JK3 coupling') then
                           write(iwrite,                               &
                             '(6x,3x,3(i2,a1,"(",i1,")"),6x,a6,f10.7)')&
                             (all_classifications%couplings(icoupling)%&
                             states(istate)%csf%subc_cfg(j)%in,CVAL(1, &
                             all_classifications%couplings(icoupling)% &
                             states(istate)%csf%subc_cfg(j)%il),       &
                             all_classifications%couplings(icoupling)% &
                             states(istate)%csf%subc_cfg(j)%iN_big,j=1,&
                             3),'coeff:',all_classifications%couplings(&
                           icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'cLSJ3 coupling') then
                        number_1 = all_classifications%couplings(    &
                           icoupling)%states(istate)%                &
                        csf%iM1(3)
                        number_2 = all_classifications%couplings(    &
                           icoupling)%states(istate)%                &
                        csf%iM2(3)
                        if(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%    &
                            subc_cfg(1)%il /=0 ) then
                   call gettermjj (2*all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%        &
                   subc_cfg(1)%il -1,all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%          &
                            iM1(1),jj_1_term,num_1)
                         else
                            jj_1_term(number_1)%j = -1
                            jj_1_term(number_1)%subshellJ = 0
                            jj_1_term(number_1)%nu  = 0
                         end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),jj_2_term,num_2)
                           write(iwrite,'(6x,3x,2(i2,a2,"(",i1,")"),&
                              2(i2,a1,"(",i1,")"),6x,a6,f10.7)')       &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),                  &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc_cfg(j)%in,          &
                           CVAL(1,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%              &
                           subc_cfg(j)%il),all_classifications%        &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(j)%iN_big, j=2,3),'coeff:',        &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%coeff
!
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj3 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(2),2,1000)
                           number_3 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),1,1000)
                           number_4 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(2),1,1000)
                           number_5 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(3),2,1000)
                           number_6 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(3),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),2,1000),jj_2_term,num_2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),1,1000),jj_3_term,num_3)
                           else
                              jj_3_term(number_3)%j = -1
                              jj_3_term(number_3)%subshellJ = 0
                              jj_3_term(number_3)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),1,1000),jj_4_term,num_4)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(3)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(3),2,1000),jj_5_term,num_5)
                           else
                              jj_5_term(number_5)%j = -1
                              jj_5_term(number_5)%subshellJ = 0
                              jj_5_term(number_5)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(3)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM1(3),1,1000),jj_6_term,num_6)
                           write(iwrite,'(9x,i2,a2,"(",i1,")",         &
                           5(i2,a2,"(",i1,")",6x),a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_3_term(number_3)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_4_term(number_4)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il,jj_5_term(number_5)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(3),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il,jj_6_term(number_6)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(3),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else
                           write(*,*) 'To many coupled shells in ',    &
                           'Coupling_evaluation.f90 ',                 &
                           'print_classification_data'
                           stop
                     end if
                  else
                     write(iwrite,*)                                   &
                        '     all_classifications%couplings()%states(',&
                        ')%csf%subc_cfg NOT associated'
                  end if
!......................................................................
                  if(associated(all_classifications%couplings(         &
                     icoupling)%states(istate)%csf%subc) .and.         &
                     associated(all_classifications%couplings(         &
                     icoupling)%states(istate)%csf%iM1) .and.          &
                     associated(all_classifications%couplings(         &
                                icoupling)%states(istate)%csf%iM2)) then
                     if(coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                             icoupling)%long_name == 'LS coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite, '(13x,i2,a1,i1)')             &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 2)  &
                                                                    then
                           write(iwrite, '(12x,7(1x,i2,a1,i1))')       &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iS+1,CVAL(2,all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              iM2(2)+1,CVAL(2,all_classifications%     &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(2))
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,3(1x,i2,a1,i1),3(1x,i2,a1))')      &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iS+1,CVAL(2,all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,3),                      &
                              (all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              iM2(j)+1,CVAL(2,all_classifications%     &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(j)), j=2,3)
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                                        long_name == 'jj1 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]" &
                                                    )')                &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'JJ coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%nosubc == 1) then
                           write(iwrite, '(13x,i2,a1,i1)')             &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 2)  &
                                                                    then
                           write(iwrite, '(12x,2(1x,i2,a1,i1),3x,"(",  &
                              a4,",",a4,")",a4)') (all_classifications%&
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),JVAL(                 &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM2(1)),JVAL(         &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM1(2)),JVAL(         &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM2(2))
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                                        long_name == 'LK coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite, '(13x,i2,a1,i1)')             &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 2)  &
                                                                    then
                           write(iwrite, '(12x,2(1x,i2,a1,i1),3x,a1,   &
                              " [",a4,"]",a4)') (all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),CVAL(2,               &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM1(2)),JVAL(         &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM2(2)),JVAL(states%  &
                              states(istate)%J)
                        end if

                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'JK coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite, '(13x,i2,a1,i1)') &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(1)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(1)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc==2) then
                           write(iwrite,                               &
                          '(12x,2(1x,i2,a1,i1),3x,a4," [",a4,"]",a4)') &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc(j)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(j)%inr, j=1,2),     &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(1)),      &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2)),      &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(2))
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LScjj coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(13x,i2,a1,i1,                &
                           2(a4,i2,1x)," [",a4,","a4," ]")')           &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr,          &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2))
                        else if(all_classifications%couplings(         &
                         icoupling)%states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(12x,2(1x,i2,a1,i1),3x,a1,    &
                              " [",a4,"]",a4)') (all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),CVAL(2,               &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM1(2)),              &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(2)),   &
                              JVAL(states%states(istate)%J)
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                                        long_name == 'jj2 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]",&
                           2(a4,1x,i1,1x,"[",a4,"]"))')                &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(jj_3_term(number_3)%subshellJ),     &
                              jj_3_term(number_3)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2)),    &
                              JVAL(jj_4_term(number_4)%subshellJ),     &
                              jj_4_term(number_4)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if

                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'LS3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,                               &
                          '(12x,i2,a1,i1," (",i2,a1,i1,1x,i2,a1,i1,")",&
                          i2,a1,1x,i2,a1)')                           &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc(j)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(j)%inr, j=1,3),     &
                           (all_classifications%couplings(             &
                           icoupling)%states(istate)%csf%iM2(j)+1,     &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%iM2(j)),j=2,3)
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'LSJ3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,                               &
                           '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,i2,a1,i1,&
                           ")",i2,a1,1x," [",a4,","a4," ]")')          &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc(j)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(j)%inr, j=1,3),     &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%iM2(2)+1,     &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%iM1(1)),      &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iJ(1)),       &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iJ(2))
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'LK3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,                               &
                          '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,i2,a1,i1, &
                          ")",i2,a1,1x,a1," [",a4," ]")')              &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc(j)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(j)%inr, j=1,3),     &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%iM2(2)+1,     &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%iM1(1)),      &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%iM1(3)),      &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(2))
                        end if
                     else if(coupling_descriptions(all_expansions%     &
                        coupling_expansions(icoupling)%icoupling)%     &
                        long_name == 'JK3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,                               &
                          '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,i2,a1,i1, &
                          ")",i2,a1,1x,a4,";",1x,a4,"[",a4" ]")')      &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc(j)%iS+1,            &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%iL),  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc(j)%inr, j=1,3),     &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%iM2(2)+1,     &
                           CVAL(2,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%iM1(2)),      &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iJ(1)),       &
                           JVAL(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iJ(2)),       &
                           JVAL(states%states(istate)%J)
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'cLSJ3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,                               &
                                 '(10x,2(1x,a4,i2)," (",i2,a1,i1,1x,i2,&
                                 a1,i1,")",i2,a1,1x," [",a4,","a4,     &
                                 " ]")')                               &
                                 JVAL(jj_1_term(number_1)%subshellJ),  &
                                 jj_1_term(number_1)%nu,               &
                                 JVAL(jj_2_term(number_2)%subshellJ),  &
                                 jj_2_term(number_2)%nu,               &
                                 (all_classifications%couplings(       &
                                 icoupling)%states(istate)%csf%subc(j)%&
                                 iS+1,CVAL(2,all_classifications%      &
                                 couplings(icoupling)%states(istate)%  &
                                 csf%subc(j)%iL),                      &
                                 all_classifications%couplings(        &
                                 icoupling)%states(istate)%csf%subc(j)%&
                                 inr, j=2,3),all_classifications%      &
                                 couplings(icoupling)%states(istate)%  &
                                 csf%iM2(2)+1,CVAL(2,                  &
                                 all_classifications%couplings(        &
                                 icoupling)%states(istate)%csf%iM1(2)),&
                                 JVAL(all_classifications%couplings(   &
                                 icoupling)%states(istate)%csf%iJ(1)), &
                                 JVAL(all_classifications%couplings(   &
                                 icoupling)%states(istate)%csf%iJ(2))
                        end if
!
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj3 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 3) then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]",&
                           4(a4,1x,i1,1x,"[",a4,"]"))')                &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(jj_3_term(number_3)%subshellJ),     &
                              jj_3_term(number_3)%nu,                  &
                              JVAL(JTHN(all_classifications%couplings( &
                              icoupling)%states(istate)%csf%iJ(2),     &
                              2,1000)),                                &
                              JVAL(jj_4_term(number_4)%subshellJ),     &
                              jj_4_term(number_4)%nu,                  &
                              JVAL(JTHN(all_classifications%couplings( &
                              icoupling)%states(istate)%csf%iJ(2),     &
                              1,1000)),                                &
                              JVAL(jj_5_term(number_5)%subshellJ),     &
                              jj_5_term(number_5)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(3)),    &
                              JVAL(jj_6_term(number_6)%subshellJ),     &
                              jj_6_term(number_6)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if
                     end if
                  else
                    write(iwrite,*) '     ',                           &
                    'all_classifications%couplings()%states()%csf%subc'&
                    ,' iM1, iM2 NOT associated'
                  end if
               end do
            else
               write(iwrite, *)                                       &
                'all_classifications%couplings()%states NOT associated'
            end if
         end do
      else
      write(iwrite, *) 'all_classifications%couplings NOT assoccated'
      end if
!
!      write(*,*) '    end subroutine print_classification_data'
!
      end subroutine print_classification_data
!
!--------------------------------------------------------------------
!
      subroutine print_classification(iwrite, icoupling)
!--------------------------------------------------------------------
!     This subroutine prints classification data
!     of coupling with serial number "icoupling" to
!     the unit "iwrite"
!--------------------------------------------------------------------
      implicit none
      integer, intent(in):: iwrite, icoupling
      integer:: istate, isubc, j
      integer:: iN_big, in, il, Li, Si, inr, L_i, S_i
      integer          :: num_1, num_2, number_1, number_2
      integer          :: num_3, num_4, number_3, number_4
      integer          :: num_5, num_6, number_5, number_6
      character(len=1) :: CVAL
      character(len=2) :: LSJVAL
      character(len=4) :: JVAL
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
      type(subshell_term), dimension(1:63)  :: jj_5_term, jj_6_term
      if(associated(all_classifications%couplings)) then
         if(all_classifications%nr_of_couplings.ge.icoupling) then
            write(iwrite,*)' '
            write(iwrite,*)      &
                     ' ---------------------------------------------- '
            write(iwrite,*)' '
            write(iwrite,*)'coupling: ', icoupling, ' (',             &
               coupling_descriptions(all_expansions%                   &
               coupling_expansions(icoupling)%icoupling)%long_name, ')'
            write(iwrite,*)' '
            if(all_classifications%couplings(icoupling)%is_one_to_one) &
                                                                    then
               write(iwrite,*)'   classification is one-to-one'
            else
               write(iwrite,*)'   classification is NOT one-to-one'
            end if
            write(iwrite,*)' '
            if(associated(all_classifications%couplings(icoupling)%    &
                                                           states)) then
               do istate=1,all_classifications%couplings(icoupling)%   &
                                                          nr_of_states,1
                  write(iwrite,*)'   ----------------------------'
                  write(iwrite,*)'   state: ', istate
                  write(iwrite,'(5x,a3,a4)')'J =', jval(states%        &
                     states(istate)%J)
                  write(iwrite,'(5x,a8,f15.8)')'Energy =',states%      &
                     states(istate)%energy
                  write(iwrite,*)' '
                  if(associated(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%subc_cfg)) then
                     if(coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj1 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(1),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(1),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM1(1),1,1000),jj_2_term,num_2)
                           write(iwrite,                               &
                           '(9x,i2,a2,"(",i1,")",                      &
                           i2,a2,"(",i1,")",6x,a6,f10.7)')             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff

                     else if(all_classifications%couplings(icoupling)% &
                                    states(istate)%csf%nosubc == 1) then
                        write(iwrite,'(1x,3x,i2,a1,"(",i1,")",6x,a6,   &
                           f10.7)') all_classifications%couplings(     &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,CVAL(1,all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il),all_classifications%couplings(          &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           iN_big,'coeff:',all_classifications%        &
                           couplings(icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LScjj coupling') then
                           number_1 = all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM1(2)
                           number_2 = all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),jj_2_term,num_2)
                           write(iwrite,                               &
                           '(9x,i2,a1,"(",i1,")",2x,             &
                           2(i2,a2,"(",i1,")"),6x,a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           CVAL(1,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il),all_classifications%couplings(          &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           iN_big,                                     &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(2)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_1_term(number_1)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(2)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_2_term(number_2)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),                  &
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj2 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),2,1000)
                           number_3 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),1,1000)
                           number_4 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),2,1000),jj_2_term,num_2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),1,1000),jj_3_term,num_3)
                           else
                              jj_3_term(number_1)%j = -1
                              jj_3_term(number_1)%subshellJ = 0
                              jj_3_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),1,1000),jj_4_term,num_4)
                           write(iwrite,                               &
                           '(9x,i2,a2,"(",i1,")",                &
                           3(i2,a2,"(",i1,")",6x),a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_3_term(number_3)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_4_term(number_4)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff

                     else if(all_classifications%couplings(icoupling)% &
                                    states(istate)%csf%nosubc == 2) then
                        write(iwrite,'(1x,3x,2(i2,a1,"(",i1,")"),6x,a6,&
                           f10.7)') (all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(j)%  &
                           in,CVAL(1,all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(j)%  &
                           il),all_classifications%couplings(          &
                           icoupling)%states(istate)%csf%subc_cfg(j)%  &
                           iN_big,j=1,2),'coeff:',all_classifications% &
                           couplings(icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LS3 coupling' .or.    &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LS coupling' .or.     &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LSJ3 coupling' .or.   &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                           icoupling)%long_name == 'LK3 coupling' .or. &
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                           icoupling)%long_name == 'JK3 coupling') then
                           write(iwrite,                               &
                             '(6x,3x,3(i2,a1,"(",i1,")"),6x,a6,f10.7)')&
                             (all_classifications%couplings(icoupling)%&
                             states(istate)%csf%subc_cfg(j)%in,CVAL(1, &
                             all_classifications%couplings(icoupling)% &
                             states(istate)%csf%subc_cfg(j)%il),       &
                             all_classifications%couplings(icoupling)% &
                             states(istate)%csf%subc_cfg(j)%iN_big,j=1,&
                             3),'coeff:',all_classifications%couplings(&
                             icoupling)%states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'cLSJ3 coupling') then
                        number_1 = all_classifications%couplings(    &
                           icoupling)%states(istate)%                &
                        csf%iM1(3)
                        number_2 = all_classifications%couplings(    &
                           icoupling)%states(istate)%                &
                        csf%iM2(3)
                        if(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%    &
                            subc_cfg(1)%il /=0 ) then
                   call gettermjj (2*all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%        &
                   subc_cfg(1)%il -1,all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%          &
                            iM1(1),jj_1_term,num_1)
                         else
                            jj_1_term(number_1)%j = -1
                            jj_1_term(number_1)%subshellJ = 0
                            jj_1_term(number_1)%nu  = 0
                         end if
                         call gettermjj(2*all_classifications%         &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),jj_2_term,num_2)
!                           write(iwrite,'(1x,i5,3x,2(i2,a2,"(",i1,")"),&
                           write(iwrite,'(9x,2(i2,a2,"(",i1,")"),     &
                              2(i2,a1,"(",i1,")"),6x,a6,f10.7)')       &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),                  &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),                 &
                           (all_classifications%couplings(icoupling)%  &
                           states(istate)%csf%subc_cfg(j)%in,          &
                           CVAL(1,all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%              &
                           subc_cfg(j)%il),all_classifications%        &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(j)%iN_big, j=2,3),'coeff:',        &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%coeff
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj3 coupling') then
                           number_1 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),2,1000)
                           number_2 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),2,1000)
                           number_3 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM1(2),1,1000)
                           number_4 =                                  &
                           JTHN(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iM2(2),1,1000)
                           number_5 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(3),2,1000)
                           number_6 =                                  &
                           JTHN(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%iM2(3),1,1000)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(1)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),2,1000),jj_1_term,num_1)
                           else
                              jj_1_term(number_1)%j = -1
                              jj_1_term(number_1)%subshellJ = 0
                              jj_1_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(1)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),2,1000),jj_2_term,num_2)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(2)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(1),1,1000),jj_3_term,num_3)
                           else
                              jj_3_term(number_1)%j = -1
                              jj_3_term(number_1)%subshellJ = 0
                              jj_3_term(number_1)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(2)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM2(1),1,1000),jj_4_term,num_4)
                           if(all_classifications%couplings(           &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il /=0 ) then
                              call gettermjj(2*all_classifications%    &
                              couplings(icoupling)%states(istate)%csf% &
                              subc_cfg(3)%il-1,                        &
                              JTHN(all_classifications%                &
                              couplings(icoupling)%states(istate)%csf% &
                              iM1(3),2,1000),jj_5_term,num_5)
                           else
                              jj_5_term(number_5)%j = -1
                              jj_5_term(number_5)%subshellJ = 0
                              jj_5_term(number_5)%nu  = 0
                           end if
                           call gettermjj(2*all_classifications%       &
                           couplings(icoupling)%states(istate)%csf%    &
                           subc_cfg(3)%il+1,                           &
                           JTHN(all_classifications%                   &
                           couplings(icoupling)%states(istate)%csf%    &
                           iM1(3),1,1000),jj_6_term,num_6)
                           write(iwrite,'(9x,i2,a2,"(",i1,")",         &
                           5(i2,a2,"(",i1,")",6x),a6,f10.7)')          &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%subc_cfg(1)%in,          &
                           LSJVAL(all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_1_term(number_1)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(1)%  &
                           il,jj_2_term(number_2)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_3_term(number_3)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(2)%  &
                           il,jj_4_term(number_4)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM2(1),1,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il,jj_5_term(number_5)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(3),2,1000),          &
                           all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           in,LSJVAL(all_classifications%couplings(    &
                           icoupling)%states(istate)%csf%subc_cfg(3)%  &
                           il,jj_6_term(number_6)%j),JTHN(             &
                           all_classifications%couplings(icoupling)%   &
                           states(istate)%csf%iM1(3),1,1000),          &
!
                           'coeff:',all_classifications%couplings(     &
                           icoupling)%states(istate)%coeff
                     else
                        write(*,*)                                     &
                  'To many coupled shells in Coupling_evaluation.f90 ',&
                        coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name
                        stop
                     end if
                  else
                     write(iwrite,*)                                   &
                       '     all_classifications%couplings()%states(', &
                       ')%csf%subc_cfg NOT associated'
                  end if
!......................................................................
                  if(associated(all_classifications%couplings(         &
                     icoupling)%states(istate)%csf%subc) .and.         &
                     associated(all_classifications%couplings(         &
                     icoupling)%states(istate)%csf%iM1) .and.          &
                     associated(all_classifications%couplings(         &
                                icoupling)%states(istate)%csf%iM2)) then
!                     if(icoupling == 1) then
                     if(coupling_descriptions(                         &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LS coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(13x,i2,a1,i1)')              &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 2)  &
                                                                    then
                           write(iwrite,'(12x,7(1x,i2,a1,i1))')        &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,2),all_classifications%couplings(&
                              icoupling)%states(istate)%csf%iM2(2)+1,  &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(2))
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,3(1x,i2,a1,i1),3(1x,i2,a1))')      &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,3),                              &
                              (all_classifications%couplings(&
                              icoupling)%states(istate)%csf%iM2(j)+1,  &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(j)),   &
                              j=2,3)
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj1 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]",&
                           a4,1x,i1,1x,"[",a4,"]")')                   &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'JJ coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(13x,i2,a1,i1)')              &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                              icoupling)%states(istate)%csf%nosubc== 2)&
                                                                    then
                           write(iwrite,'(12x,2(1x,i2,a1,i1),3x,"(",a4,&
                              ",",a4,")",a4)') (all_classifications%   &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),                      &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(1)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM1(2)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(2))
                        end if
!                     else if(icoupling == 3) then
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LK coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(13x,i2,a1,i1)')              &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                         icoupling)%states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(12x,2(1x,i2,a1,i1),3x,a1,    &
                              " [",a4,"]",a4)') (all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),CVAL(2,               &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%iM1(2)),              &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(2)),   &
                              JVAL(states%states(istate)%J)
                        end if

!                     else if(icoupling == 4) then
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'JK coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 1) then
                           write(iwrite,'(13x,i2,a1,i1)')              &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr
                        else if(all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%nosubc == 2)  &
                                                                    then
                           write(iwrite,'(12x,2(1x,i2,a1,i1),3x,a4,    &
                              " [",a4,"]",a4)') (all_classifications%  &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iS+1,CVAL(2,all_classifications% &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%iL),all_classifications%         &
                              couplings(icoupling)%states(istate)%csf% &
                              subc(j)%inr,j=1,2),                      &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(1)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM1(2)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(2))
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LScjj coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(13x,i2,a1,i1,                &
                           2(a4,i2,1x)," [",a4,","a4," ]")')           &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iS+1,CVAL(2,  &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%iL),          &
                              all_classifications%couplings(icoupling)%&
                              states(istate)%csf%subc(1)%inr,          &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2))
                        end if
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'jj2 coupling') then
                        if(all_classifications%couplings(icoupling)%   &
                                    states(istate)%csf%nosubc == 2) then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]",&
                           2(a4,1x,i1,1x,"[",a4,"]"))')                &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(jj_3_term(number_3)%subshellJ),     &
                              jj_3_term(number_3)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2)),    &
                              JVAL(jj_4_term(number_4)%subshellJ),     &
                              jj_4_term(number_4)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if

!                     else if(icoupling == 5) then
                      else if(coupling_descriptions(                   &
                         all_expansions%coupling_expansions(icoupling)%&
                         icoupling)%long_name == 'LS3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,i2,a1,i1," (",i2,a1,i1,1x,i2,a1,   &
                              i1,")",i2,a1,1x,i2,a1)')                 &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,3),                              &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%iM2(j)+1,  &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(j)),   &
                              j=2,3)
                        end if
!                     else if(icoupling == 6) then
                      else if(coupling_descriptions(                   &
                         all_expansions%coupling_expansions(icoupling)%&
                         icoupling)%long_name == 'LSJ3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,      &
                              i2,a1,i1,")",i2,a1,1x," [",a4,","a4,     &
                              " ]")')                                  &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,3),                              &
                              all_classifications%couplings(           &
                              icoupling)%states(istate)%csf%iM2(2)+1,  &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(2)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2))
                        end if
!                     else if(icoupling == 7) then
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'LK3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,i2,a1,&
                              i1,")",i2,a1,1x,a1," [",a4," ]")')       &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,3),                              &
                              all_classifications%couplings(           &
                              icoupling)%states(istate)%csf%iM2(2)+1,  &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(2)),   &
                              CVAL(2,all_classifications%couplings(    &
                              icoupling)%states(istate)%csf%iM1(3)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM2(3))
                        end if
!                     else if(icoupling == 8) then
                     else if(coupling_descriptions(                    &
                        all_expansions%coupling_expansions(icoupling)% &
                        icoupling)%long_name == 'JK3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x,      &
                              i2,a1,i1,")",i2,a1,1x,a4,";",1x,a4," [", &
                              a4," ]")')                               &
                              (all_classifications%couplings(          &
                              icoupling)%states(istate)%csf%subc(j)%iS+&
                              1,CVAL(2,all_classifications%couplings(  &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              iL),all_classifications%couplings(       &
                              icoupling)%states(istate)%csf%subc(j)%   &
                              inr,j=1,3),                              &
                              all_classifications%couplings(           &
                              icoupling)%states(istate)%csf%iM2(2)+1,  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iM1(2)),   &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(2)),    &
                              Jval(states%states(istate)%J)
                        end if
!                     else if(icoupling == 9) then
                     else if(coupling_descriptions(                   &
                         all_expansions%coupling_expansions(icoupling)%&
                         icoupling)%long_name == 'cLSJ3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,                               &
                              '(10x,2(1x,a4,i2)," (",i2,a1,i1,1x,i2,&
                              a1,i1,")",i2,a1,1x," [",a4,","a4,     &
                              " ]")')                               &
                           JVAL(jj_1_term(number_1)%subshellJ),  &
                           jj_1_term(number_1)%nu,               &
                           JVAL(jj_2_term(number_2)%subshellJ),  &
                           jj_2_term(number_2)%nu,               &
                           (all_classifications%couplings(       &
                           icoupling)%states(istate)%csf%subc(j)%&
                           iS+1,CVAL(2,all_classifications%      &
                           couplings(icoupling)%states(istate)%  &
                           csf%subc(j)%iL),                      &
                           all_classifications%couplings(        &
                           icoupling)%states(istate)%csf%subc(j)%&
                           inr, j=2,3),all_classifications%      &
                           couplings(icoupling)%states(istate)%  &
                           csf%iM2(2)+1,CVAL(2,                  &
                           all_classifications%couplings(        &
                           icoupling)%states(istate)%csf%iM1(2)),&
                           JVAL(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iJ(1)), &
                           JVAL(all_classifications%couplings(   &
                           icoupling)%states(istate)%csf%iJ(2))
                        end if
!
                     else if(coupling_descriptions(                    &
                         all_expansions%coupling_expansions(icoupling)%&
                         icoupling)%long_name == 'jj3 coupling') then
                        if(all_classifications%couplings(              &
                           icoupling)%states(istate)%csf%nosubc == 3)  &
                                                                    then
                           write(iwrite,'(9x,2(a4,1x,i1,1x),"[",a4,"]",&
                           4(a4,1x,i1,1x,"[",a4,"]"))')                &
                              JVAL(jj_1_term(number_1)%subshellJ),     &
                              jj_1_term(number_1)%nu,                  &
                              JVAL(jj_2_term(number_2)%subshellJ),     &
                              jj_2_term(number_2)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(1)),    &
                              JVAL(jj_3_term(number_3)%subshellJ),     &
                              jj_3_term(number_3)%nu,                  &
                              JVAL(JTHN(all_classifications%couplings( &
                              icoupling)%states(istate)%csf%iJ(2),     &
                              2,1000)),                                &
                              JVAL(jj_4_term(number_4)%subshellJ),     &
                              jj_4_term(number_4)%nu,                  &
                              JVAL(JTHN(all_classifications%couplings( &
                              icoupling)%states(istate)%csf%iJ(2),     &
                              1,1000)),                                &
                              JVAL(jj_5_term(number_5)%subshellJ),     &
                              jj_5_term(number_5)%nu,                  &
                              JVAL(all_classifications%couplings(      &
                              icoupling)%states(istate)%csf%iJ(3)),    &
                              JVAL(jj_6_term(number_6)%subshellJ),     &
                              jj_6_term(number_6)%nu,                  &
                              JVAL(states%states(istate)%J)
                        end if
                     else
                        do isubc=1,all_classifications%couplings(      &
                           icoupling)%states(istate)%csf%nosubc,1
                           inr = all_classifications%couplings(        &
                           icoupling)%states(istate)%csf%subc(isubc)%inr
                           Li = all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%subc(isubc)%iL
                           Si = all_classifications%couplings(         &
                           icoupling)%states(istate)%csf%subc(isubc)%iS
                           S_i = all_classifications%couplings(        &
                           icoupling)%states(istate)%csf%iM2(isubc)
                           L_i = all_classifications%couplings(        &
                           icoupling)%states(istate)%csf%iM1(isubc)
                           write(iwrite,'(9x,i2,3x,i2,3x,i2,3x,i2,3x,  &
                              i2,3x,i2)')isubc,inr,Li,Si,L_i,S_i
                        end do
                     end if
                     write(iwrite,*)' '
                  else
                     write(iwrite,*)  &
             '     all_classifications%couplings()%states()%csf%subc,',&
                    ' iM1, iM2 NOT associated'
                  end if
               end do
            else
               write(iwrite, *)       &
               'all_classifications%couplings()%states NOT associated'
            end if
         else
            write(*,*) 'ERROR in subroutine print_classification: ',   &
               'all_classifications%nr_of_couplings.ge.icoupling; ',   &
               'classification data will not be printed !'
         end if
      else
         write(iwrite, *) 'all_classifications%couplings NOT assoccated'
      end if
      end subroutine print_classification
!
!***********************************************************************
!                                                                      *
      SUBROUTINE getchLS(shell_nn,shell_ll,shell_N,shell_L,shell_S,    &
      shell_w,string_shell)
!                                                                      *
!     The spectroscopic notation of a shell in LS coupling is returned.*
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                 last update: January 2020   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: shell_nn,shell_ll,shell_N,shell_L
      integer, intent(in)           :: shell_S, shell_w
      character(len=9), intent(out) :: string_shell
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: string_lenth, number, j, jnumber, nu
!GG      character(len=1) :: string_N,string_S,string_v
      character(len=1) :: string_S,string_v
      character(len=2) :: string_N
      character(len=2) :: string_nn
      character(len=4) :: string_CNUM
      type(subshell_term_LS), dimension(120) :: LS
!-----------------------------------------------
      call convrt(shell_nn,string_CNUM,string_lenth)
      string_nn = string_CNUM(1:string_lenth)
      call convrt_double(2*(1+shell_S),string_CNUM,string_lenth)
      string_S = string_CNUM(1:string_lenth)
      if(shell_N /=0 .or. shell_N /=1 .or. shell_N /= 4*shell_ll+2)then
         call gettermLS(shell_ll,shell_N,LS,number)
         jnumber = 0
         do j =1,number
            if(shell_L==LS(j)%LL .and. shell_S==LS(j)%S)               &
                                                      jnumber=jnumber+1
         end do
         if (jnumber /=1 .and. shell_ll < 3) then
            nu =                                                       &
            Seniority(shell_ll,shell_N,shell_L,shell_S,shell_w)
            call convrt_double(2*nu,string_CNUM,string_lenth)
            string_v = string_CNUM(1:string_lenth)
         else if (jnumber /=1) then
            call convrt_double(2*shell_w,string_CNUM,string_lenth)
            string_v = string_CNUM(1:string_lenth)
         end if
      end if
      if (shell_N==1) then
         string_shell = trim(string_nn)//ll_string(shell_ll)//'_'
      else if (shell_N==4*shell_ll+2) then
         call convrt(shell_N,string_CNUM,string_lenth)
         string_N = trim(string_CNUM(1:string_lenth))
         string_shell =                                               &
             trim(string_nn)//ll_string(shell_ll)//trim(string_N)//'_'
!GG                    trim(string_nn)//ll_string(shell_ll)//string_N//'_'
      else
         call convrt(shell_N,string_CNUM,string_lenth)
         string_N = string_CNUM(1:string_lenth)
         if (jnumber == 1) then
            string_shell =                                             &
               trim(string_nn)//ll_string(shell_ll)//trim(string_N)//  &
!GG                     trim(string_nn)//ll_string(shell_ll)//string_N//  &
                    '('//string_S//L_string(shell_L/2)//')'
         else
            string_shell =                                             &
               trim(string_nn)//ll_string(shell_ll)//trim(string_N)//  &
!GG                     trim(string_nn)//ll_string(shell_ll)//string_N//  &
                    '('//string_S//L_string(shell_L/2)//string_v//')'
         end if
      end if
!
      END SUBROUTINE getchLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE getchjj(shell_nn,shell_ll,shell_jj,shell_N,shell_J,   &
      shell_w,string_shell)
!                                                                      *
!     The spectroscopic notation of a shell in LS coupling is returned.*
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)          :: shell_nn,shell_ll,shell_jj,shell_N
      integer, intent(in)           :: shell_J, shell_w
      character(len=9), intent(out) :: string_shell
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: string_lenth, number, j, jnumber, nu
      character(len=1) :: string_N,string_J,string_v
      character(len=2) :: LSJVAL, string_jj
      character(len=4) :: string_CNUM
      character(len=2) :: string_nn
      type(subshell_term), dimension(63) :: jj
!-----------------------------------------------
      call convrt(shell_nn,string_CNUM,string_lenth)
      string_nn = string_CNUM(1:string_lenth)
      string_jj = lsjval(shell_ll,shell_jj)
!GG NK pradzia
      call convrt_double(shell_J,string_CNUM,string_lenth)
      string_J = string_CNUM(1:string_lenth)
      if(shell_N /=0 .or. shell_N /=1 .or. shell_N /= shell_jj+1)then
         call gettermjj(shell_jj,shell_N,jj,number)
         jnumber = 0
         do j =1,number
            if(shell_J == jj(j)%subshellJ) jnumber=jnumber+1
         end do
!NK         if (jnumber /=1 .and. shell_jj < 9) then
!NK            nu =                                                       &
!NK            Seniority_jj(shell_jj,shell_N,shell_J,shell_w)
!NK           call convrt_double(2*nu,string_CNUM,string_lenth)
!NK            string_v = string_CNUM(1:string_lenth)
!NK         else if (jnumber /=1) then
         if (jnumber /=1) then
            call convrt_double(2*shell_w,string_CNUM,string_lenth)
            string_v = string_CNUM(1:string_lenth)
         end if
      end if
!GG NK pabaiga
      if (shell_N==1) then
         string_shell = trim(string_nn)//string_jj//'_'
      else if (shell_N==0) then
         string_shell = ''
      else if (shell_N== shell_jj+1) then
         call convrt(shell_N,string_CNUM,string_lenth)
         string_N = string_CNUM(1:string_lenth)
         string_shell = trim(string_nn)//string_jj//string_N//'_'
      else
         call convrt(shell_N,string_CNUM,string_lenth)
         string_N = string_CNUM(1:string_lenth)
         string_shell = trim(string_nn)//string_jj//string_N
!GG NK pradzia
         if (jnumber /= 1) then
            string_shell = trim(string_nn)//string_jj//string_N//  &
                    '{'//string_v//'}'
!                    '<'//string_J//'>'//'('//string_v//')'
         else
!            string_shell = trim(string_nn)//string_jj//string_N//  &
!                    '<'//string_J//'>'
         end if
!GG NK pabaiga
      end if
!
      END SUBROUTINE getchjj
!
!***********************************************************************
!                                                                      *
      SUBROUTINE getchXLS(shell_L,shell_S,string_shell)
!                                                                      *
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2015   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: shell_L,shell_S
      character(len=2), intent(out) :: string_shell
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: string_lenth
      character(len=1) :: string_S
      character(len=4) :: string_CNUM
!-----------------------------------------------
      call convrt_double(2*(1+shell_S),string_CNUM,string_lenth)
      string_S = string_CNUM(1:string_lenth)
      string_shell = string_S//L_string(shell_L/2)
!
      END SUBROUTINE getchXLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE CONVRT_DOUBLE(INTNUM, CNUM, LENTH)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by G. Gaigalas,                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:53   2/14/04
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,   INTENT(IN)    :: INTNUM
      INTEGER,   INTENT(OUT)   :: LENTH
      CHARACTER, INTENT(INOUT) :: CNUM*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                           :: INTNUMGG
      CHARACTER(LEN=6)                  :: FORM
      CHARACTER(LEN=2), DIMENSION(0:10) :: C1020
      CHARACTER, DIMENSION(9)           :: C19
!
      DATA C19 /'1','2','3','4','5','6','7','8','9'/
      DATA C1020/'10','11','12','13','14','15','16','17','18','19','20'/
!-----------------------------------------------
      IF(mod(INTNUM,2) == 0) THEN
         INTNUMGG = INTNUM/2
      ELSE
         INTNUMGG = INTNUM
      ENDIF
!
      IF (INTNUMGG < 0) THEN
         LENTH = LOG10(DBLE((-INTNUMGG))) + 2
      ELSE IF (INTNUMGG == 0) THEN
         LENTH = 1
      ELSE
         LENTH = LOG10(DBLE(INTNUMGG)) + 1
      ENDIF
!
!   Ensure that the length of CNUM as dimensioned is adequate;
!   stop with an error message if it isn't
!
      IF (LENTH > LEN(CNUM)) THEN
         WRITE (6, *) 'CONVRT_DOUBLE: Length of CNUM inadeuate.'
         STOP
      ELSE
         IF (LENTH <= 9) THEN
            FORM = '(1I'//C19(LENTH)//')'
            WRITE (CNUM(1:LENTH), FORM(1:5)) INTNUMGG
         ELSE
            FORM = '(1I'//C1020(LENTH-10)//')'
            WRITE (CNUM(1:LENTH), FORM(1:6)) INTNUMGG
         ENDIF
         IF(mod(INTNUM,2) /= 0) THEN
            IF (LENTH+2 > LEN(CNUM)) THEN
               WRITE (6, *) 'CONVRT_DOUBLE: Length of CNUM inadeuate.'
               STOP
            ELSE
               CNUM(1:LENTH+2) = CNUM(1:LENTH)//'/2'
               LENTH = LENTH + 2
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CONVRT_DOUBLE
!
!***********************************************************************
!                                                                      *
      SUBROUTINE CONVRT(INTNUM, CNUM, LENTH)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 22 Sep 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:53   2/14/04
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,   INTENT(IN)    :: INTNUM
      INTEGER,   INTENT(OUT)   :: LENTH
      CHARACTER, INTENT(INOUT) :: CNUM*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=6)                  :: FORM
      CHARACTER(LEN=2), DIMENSION(0:10) :: C1020
      CHARACTER, DIMENSION(9)           :: C19
!
      DATA C19 /'1','2','3','4','5','6','7','8','9'/
      DATA C1020/'10','11','12','13','14','15','16','17','18','19','20'/
!-----------------------------------------------
      IF (INTNUM < 0) THEN
         LENTH = LOG10(DBLE((-INTNUM))) + 2
      ELSE IF (INTNUM == 0) THEN
         LENTH = 1
      ELSE
         LENTH = LOG10(DBLE(INTNUM)) + 1
      ENDIF
!
!   Ensure that the length of CNUM as dimensioned is adequate;
!   stop with an error message if it isn't
!
      IF (LENTH > LEN(CNUM)) THEN
         WRITE (6, *) 'CONVRT: Length of CNUM inadeuate.'
         STOP
      ELSE
         IF (LENTH <= 9) THEN
            FORM = '(1I'//C19(LENTH)//')'
            WRITE (CNUM(1:LENTH), FORM(1:5)) INTNUM
         ELSE
            FORM = '(1I'//C1020(LENTH-10)//')'
            WRITE (CNUM(1:LENTH), FORM(1:6)) INTNUM
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE CONVRT
!
!***********************************************************************
!                                                                      *
      function Seniority (ll,N,L,S,nu)          result(irez)
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  Oct 2015   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
      integer, intent(in) :: ll, N, L, S, nu
      integer :: irez
      irez = nu
      if(NOTATION == 1) then
         if (ll == 2) then
            if (N == 3 .or. N == 7) then
!               if(L == 2) then
               if(L == 4) then
                  if (nu == 1) irez = 1
                  if (nu == 3) irez = 2
               end if
            else if (N == 4 .or. N == 6) then
               if (L == 0) then
                  if (nu == 0) irez = 1
                  if (nu == 4) irez = 2
               else if (L == 1) then
                  irez = nu/2
               else if (L == 2 .and. S == 0) then
                  irez = nu/2
               else if (L == 3 .and. S == 2) then
                  irez = nu/2
               end if
            else if (N == 5) then
!               if(L == 2 .and. S == 1) then
               if(L == 4 .and. S == 1) then
                  if (nu == 1) irez = 1
                  if (nu == 3) irez = 2
                  if (nu == 5) irez = 3
!               else if(L == 3 .and. S == 1) then
               else if(L == 6 .and. S == 1) then
                  if (nu == 3) irez = 1
                  if (nu == 5) irez = 2
!               else if(L == 4 .and. S == 1) then
               else if(L == 8 .and. S == 1) then
                  if (nu == 3) irez = 1
                  if (nu == 5) irez = 2
               end if
            end if
         end if
      end if
      end function Seniority
!
!***********************************************************************
!                                                                      *
      function Seniority_JJ (jj,N,J,nu)          result(irez)
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  Oct 2015   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
      integer, intent(in) :: jj, N, J, nu
      integer :: irez
      irez = nu
      print*,"jj",jj,"i N=",N," J=",J," nu=",nu
      if(NOTATION == 1) then
         if (jj == 7) then
            if (N == 4) then
               if(J == 2) then
                  if (nu == 2) irez = 1
                  if (nu == 4) irez = 2
               else if(J == 4) then
                  if (nu == 2) irez = 1
                  if (nu == 4) irez = 2
               end if
            end if
         else if (jj == 9) then
            if (N == 3 .or. N == 7) then
               if (J == 9) then
                  if (nu == 1) irez = 1
                  if (nu == 3) irez = 2
               end if
            else if (N == 4 .or. N == 6) then
               if(J == 0) then
                  if (nu == 0) irez = 1
                  if (nu == 4) irez = 2
               else if(J == 2) then
                  if (nu == 2) irez = 1
                  if (nu == 4) irez = 2
               else if(J == 4) then
                  if (nu == 2) irez = 1
                  if (nu == 4) irez = 2
! nepabaigta
               end if
            end if
         end if
      end if
      end function Seniority_JJ
!
      end module Coupling_evaluation
