      module Coupling_main
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
      use Coupling_getdata_mchf
      use Coupling_transform_LSJJ
      use Coupling_transform_LSLK
      use Coupling_transform_JJJK
      use Coupling_transform_LSLS3
      use Coupling_transform_LSLSJ3
      use Coupling_transform_LSLK3
      use Coupling_transform_LSJK3
      use Coupling_transform_LScLSJ3
      use Coupling_transform_LSLScjj
      use Coupling_transform_LSjj1
      use Coupling_transform_LSjj2
      use Coupling_transform_LSjj3
      use Coupling_evaluation
      use Coupling_inside_shell
!
      public  :: main
      public  :: evaluate_couplings
      private :: find_nr_of_couplings
      public  :: getdata
      public  :: printas
      public  :: print_DataBasis
      public  :: iseiti
      public  :: store_R_values
      public  :: store_P_values
!
      contains
!--------------------------------------
!
      subroutine main(idata_source, couplings, print_level,thresh)
!--------------------------------------------------------------------
! This is managerial subroutine
! controlling the actions of the data input,
! procession and output
!--------------------------------------------------------------------
!
      integer, intent(in):: idata_source
      integer, intent(in):: print_level
      character(len=13), intent(in) :: couplings
      real(kind=dp), intent(in) :: thresh
      integer:: iLS_coupling_nr, iJJ_coupling_nr
      integer:: nr_of_j, nr_of_j_count
!
      integer::icoupling
!
!      write(*,*) ' subroutine main'
      open(2000+iread_csf, status='scratch')
      iLS_coupling_nr = 0
      iJJ_coupling_nr = 0
      nr_of_j_count = 0
      nr_of_j = 1
      call fill_up_coupling_descriptions
      DO
      nr_of_j_count = nr_of_j_count + 1
      if (nr_of_j_count > nr_of_j) exit
      if(evaluate_couplings(couplings)) then
         if(find_nr_of_couplings(couplings).gt.0) then
            nr_of_requested_couplings = find_nr_of_couplings(couplings)
            all_expansions%nr_of_coupling_expansions = &
                                              nr_of_requested_couplings
            allocate(all_expansions%coupling_expansions(               &
                              all_expansions%nr_of_coupling_expansions))
!
            icoupling=0
            inrcoupling=1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: ',               &
                 'icoupling.gt.all_expansions%nr_of_coupling_expansions'
                 stop
              end if
!GG              call getdata(idata_source, icoupling, print_level)
              IF(nr_of_j_count == 1) call getdataGG(idata_source,      &
                                         icoupling,print_level,nr_of_j)
              IF(nr_of_j_count > 1) call getdata(idata_source,         &
                           icoupling,print_level,nr_of_j,nr_of_j_count)

              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              iLS_coupling_nr = icoupling
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                  all_expansions%nr_of_coupling_expansions) then
                  write(*,*) 'STOP at subroutine main: icoupling.gt.', &
                  'all_expansions%nr_of_coupling_expansions'
                  stop
              end if
              write(*,*) '  Performing the transformation LS -> JJ'
              call main_LSJJ(icoupling, 1, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling &
                 = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
              iJJ_coupling_nr = icoupling
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                 'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> LK'
              call main_LSLK(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation JJ -> JK'
              call main_JJJK(icoupling, iJJ_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> LS3'
              call main_LSLS3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> LSJ3'
              call main_LSLSJ3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> LK3'
              call main_LSLK3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> JK3'
              call main_LSJK3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> cLSJ3'
              call main_LScLSJ3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> LScjj'
              call main_LSLScjj(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> jj'
              call main_LSjj1(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> jj'
              call main_LSjj2(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
            if(couplings(inrcoupling:inrcoupling).eq.'y') then
              icoupling=icoupling+1
              if(icoupling .gt.                                        &
                 all_expansions%nr_of_coupling_expansions) then
                 write(*,*) 'STOP at subroutine main: icoupling.gt.',  &
                      'all_expansions%nr_of_coupling_expansions'
                 stop
              end if
              write(*,*) '  Performing the transformation LS -> jj'
              call main_LSjj3(icoupling, iLS_coupling_nr, print_level)
              all_expansions%coupling_expansions(icoupling)%icoupling  &
              = inrcoupling
              write(*,*) '             The transformation completed'
              write(*,*) ''
            end if
            inrcoupling = inrcoupling +1
!
! here to add all the aditional coupling types ...
!
! here to add all the aditional coupling types ...
!
! here to add all the aditional coupling types ...
!
            write(*,*) '  All transformations completed'
            write(*,*) ' '
            if(print_level.gt.0) call printas
            if(print_level.ge.0) call print_DataBasis(nr_of_j_count,thresh)
!
!classification and evaluation
            call main_evaluation(nr_of_j_count,nr_of_j,                &
            iwrite_classifications,iwrite_classifications_evaluations, &
                                         iwrite_suggestions,print_level)
         end if
      else
         stop 'STOP at subroutine main: not valid list of couplings'
      end if
      END DO
      close(2000+iread_csf)
!
!      write(*,*) ' end subroutine main'
      end subroutine main
!
!--------------------------------------
!
      function evaluate_couplings(couplings)    result(rez)
!--------------------------------------------------------------------
! This function checks the consistency of the
! input string "couplings" i.e. checks whether
! it contains the "y" and "n" characters only
!
!--------------------------------------------------------------------
!
      character(len=13), intent(in) :: couplings
      logical::rez
      integer::i
      rez=.true.
      do i=1, NR_OF_AVIALABLE_COUPLINGS, 1
         if(couplings(i:i).ne.'y'.and.couplings(i:i).ne.'n') then
            rez=.false.
            exit
         end if
      end do
      end function evaluate_couplings
!
!--------------------------------------
!
      function find_nr_of_couplings(couplings)    result(irez)
!--------------------------------------------------------------------
! This function finds the number of
! requested couplings i.e. the number of
! "y" characters in the input string "couplings"
!
!--------------------------------------------------------------------
!
      character(len=13), intent(in) :: couplings
      integer::irez
      integer::i
      irez=0
      do i=1, NR_OF_AVIALABLE_COUPLINGS, 1
         if(couplings(i:i).eq.'y') irez=irez+1
      end do
      end function find_nr_of_couplings
!
!--------------------------------------
!
      subroutine getdataGG(data_source,icoupling_nr_in_expansions,    &
                                                  print_level,nr_of_j)
!--------------------------------------------------------------------
! This processes the data input (will have sence
! when there will be more than one data input interface)
!
!--------------------------------------------------------------------
      integer,intent(in) :: data_source, icoupling_nr_in_expansions
      integer,intent(in) :: print_level
      integer,intent(out) :: nr_of_j
!      write(*,*) '  subroutine getdata'
      select case (data_source)
      case(1)
!         write(*,*)'   mchf data'
!GG         call get_mchf_data(icoupling_nr_in_expansions, print_level)
         call get_mchf_data(1, icoupling_nr_in_expansions,     &
                               print_level, nr_of_j)
      case default
         write(*,*)'Unknown data source', data_source
         write(*,*)'STOP at subroutine getdata, module Coupling_main:',&
                   ' Unknown data source'
         stop
      end select
!      write(*,*) '  end subroutine getdata'
      end subroutine getdataGG
!
!--------------------------------------
!
      subroutine getdata(data_source,icoupling_nr_in_expansions,      &
                                             print_level,nr_of_j,icase)
!--------------------------------------------------------------------
! This processes the data input (will have sence
! when there will be more than one data input interface)
!
!--------------------------------------------------------------------
      integer,intent(in) :: data_source, icoupling_nr_in_expansions
      integer,intent(in) :: print_level, icase
      integer,intent(inout) :: nr_of_j
!      write(*,*) '  subroutine getdata'
      select case (data_source)
      case(1)
!         write(*,*)'   mchf data'
!GG         call get_mchf_data(icoupling_nr_in_expansions, print_level)
         call get_mchf_data(icase, icoupling_nr_in_expansions,     &
                               print_level,nr_of_j)
      case default
         write(*,*)'Unknown data source', data_source
         write(*,*)'STOP at subroutine getdata, module Coupling_main:',&
                   ' Unknown data source'
         stop
      end select
!      write(*,*) '  end subroutine getdata'
      end subroutine getdata

!
      subroutine iseiti
!--------------------------------------------------------------------
! This subroutine deallocates the general data
! arrays
!
!--------------------------------------------------------------------
      implicit none
      integer::error
      integer::icoupling, istate, icsf
!write(*,*) ''
!      write(*,*) ' subroutine iseiti'
!
!deallocate states
!
      if(associated(states%states)) deallocate(states%states,STAT=error)
!
!deallocate expansions
!
      if(associated(all_expansions%coupling_expansions)) then
         do icoupling=1, all_expansions%nr_of_coupling_expansions
            if(associated(all_expansions%coupling_expansions(          &
                                           icoupling)%expansions)) then
               do istate=1,all_expansions%coupling_expansions(         &
                                             icoupling)%nr_of_expansions
                  if(associated(all_expansions%coupling_expansions(    &
                     icoupling)%expansions(istate)%coeffs))            &
                     deallocate(all_expansions%coupling_expansions(    &
                        icoupling)%expansions(istate)%coeffs,STAT=error)
!set up deallocation of arrays inside csfs
                  if(associated(all_expansions%coupling_expansions(    &
                              icoupling)%expansions(istate)%csfs))  then
                     do icsf=1,all_expansions%coupling_expansions(     &
                                   icoupling)%expansions(istate)%size,1
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%subc))                   &
                           deallocate(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%subc)
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%subc_cfg))               &
                           deallocate(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%subc_cfg)
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iM1))                    &
                           deallocate(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iM1)
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iM2))                    &
                           deallocate(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iM2)
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iJ))                     &
                           deallocate(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                           istate)%csfs(icsf)%iJ)
                     end do
                     deallocate(all_expansions%coupling_expansions(    &
                          icoupling)%expansions(istate)%csfs,STAT=error)
                  end if
               end do
               deallocate(all_expansions%coupling_expansions(          &
                  icoupling)%expansions, STAT = error)
            end if
         end do
      end if
!
!-----  deallocate classification data    --------------------
!
      if(associated(all_classifications%couplings)) then
         do icoupling=1, all_classifications%nr_of_couplings
            if(associated(all_classifications%couplings(icoupling)%    &
               states)) deallocate(all_classifications%couplings(      &
               icoupling)%states, STAT=error)
         end do !icoupling
         deallocate(all_classifications%couplings, STAT=error)
      end if !associated(all_classifications%couplings)
!
!-----  deallocate classification evaluation data    --------------------
!
      if(associated(all_evaluations%R)) then
         do icoupling=1, all_evaluations%nr_of_couplings
!GG            if(associated(all_evaluations%R(icoupling)%RJ_values))     &
!GG           deallocate(all_evaluations%R(icoupling)%RJ_values,STAT=error)
!GG            if(associated(all_evaluations%R(icoupling)%J)) &
!GG                  deallocate(all_evaluations%R(icoupling)%J,STAT=error)
         end do !icoupling
         deallocate(all_evaluations%R, STAT=error)
      end if !associated(all_evaluations%R)
!
!      write(*,*) ' end subroutine iseiti'
!
      end subroutine iseiti
!--------------------------------------
!
      subroutine printas
!--------------------------------------------------------------------
! This subroutine prints the expansions to the file
!
!--------------------------------------------------------------------
      implicit none
      integer:: icoupling, istate, icsf, isubc, j
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
      write(iwrite_expansions,*)'-------------------------------------'
      write(iwrite_expansions,*)'Info about states'
      write(iwrite_expansions,*)'nr_of_states:', states%nr_of_states
      write(iwrite_expansions,*)' '
!      write(*,*) ' subroutine printas'
      if(associated(all_expansions%coupling_expansions)) then
         do icoupling=1, all_expansions%nr_of_coupling_expansions, 1
            write(iwrite_expansions,*) &
               '*************************************'
            write(iwrite_expansions,*)                                 &
               coupling_descriptions(all_expansions%                   &
                      coupling_expansions(icoupling)%icoupling)%       &
                      long_name
            write(iwrite_expansions,*)' '
            if(associated(all_expansions%coupling_expansions(          &
               icoupling)%expansions)) then
               if(all_expansions%coupling_expansions(icoupling)%       &
                  nr_of_expansions.ne.states%nr_of_states) then
                  write(*,*) 'stop at subroutine printas: ',           &
                     'all_expansions%coupling_expansions(',            &
                     ')%nr_of_expansions.ne.states%nr_of_states'
                  stop
               end if
               do istate=1, states%nr_of_states
                  write(iwrite_expansions,*)  &
                      '----------------------------'
                  write(iwrite_expansions,*)'state Nr.',istate
                  write(iwrite_expansions,'(3x,a3,a4,a9,3x,f15.8)')&
                  'J = ', JVAL(states%states(istate)%J),           &
                  ' Energy =',states%states(istate)%energy
                  write(iwrite_expansions,'(3x,a28,i7)')           &
                  'expansion size: ',                              &
                  all_expansions%coupling_expansions(icoupling)%   &
                            expansions(istate)%size
                  write(iwrite_expansions,*)''
                  if(associated(all_expansions%coupling_expansions(    &
                               icoupling)%expansions(istate)%csfs)) then
                     do icsf=1,all_expansions%coupling_expansions(     &
                                    icoupling)%expansions(istate)%size,1
                        if(icsf ==1)write(iwrite_expansions,*)'csf Nr.'
                        if(associated(all_expansions%                  &
                           coupling_expansions(icoupling)%expansions(  &
                                      istate)%csfs(icsf)%subc_cfg)) then
                           if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                              'jj1 coupling' .and.                     &
                              all_expansions%coupling_expansions(      &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 1) then
                                number_1 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(1),  &
                                2,1000)
                                number_2 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(1),  &
                                1,1000)
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),2,1000),jj_1_term,num_1)
                               else
                                  jj_1_term(number_1)%j = -1
                                  jj_1_term(number_1)%subshellJ = 0
                                  jj_1_term(number_1)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),1,1000),jj_2_term,num_2)
!
                                write(iwrite_expansions,               &
                                '(1x,i5,3x,i2,a2,"(",i1,")",           &
                                i2,a2,"(",i1,")",6x,a6,f10.7)')icsf,&
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_1_term(number_1)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),2,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_2_term(number_2)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),1,1000),                        &
!
                                'coeff:',all_expansions%               &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%coeffs(icsf)
                           else if(all_expansions%coupling_expansions( &
                              icoupling)%expansions(istate)%csfs(icsf)%&
                                                       nosubc == 1) then
                              write(iwrite_expansions,                 &
                             '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)')&
                               icsf,all_expansions%coupling_expansions(&
                               icoupling)%expansions(istate)%csfs(     &
                               icsf)%subc_cfg(1)%in,CVAL(1,            &
                               all_expansions%coupling_expansions(     &
                               icoupling)%expansions(istate)%csfs(     &
                               icsf)%subc_cfg(1)%il),all_expansions%   &
                               coupling_expansions(icoupling)%         &
                               expansions(istate)%csfs(icsf)%subc_cfg( &
                               1)%iN_big,'coeff:',all_expansions%      &
                               coupling_expansions(icoupling)%         &
                               expansions(istate)%coeffs(icsf)
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                              'LScjj coupling' .and.                   &
                              all_expansions%coupling_expansions(      &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 2) then
                                number_1 = all_expansions%             &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(2)
                                number_2 = all_expansions%             &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2)
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il -1, all_expansions%   &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),jj_1_term,num_1)
                               else
                                  jj_1_term(number_1)%j = -1
                                  jj_1_term(number_1)%subshellJ = 0
                                  jj_1_term(number_1)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il+1,all_expansions%       &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(1),  &
                                jj_2_term,num_2)
                                write(iwrite_expansions,               &
                                '(1x,i5,3x,i2,a1,"(",i1,")",2x,        &
                                2(i2,a2,"(",i1,")"),6x,a6,f10.7)')icsf,&
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                CVAL(1,all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il),all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%iN_big,                    &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_1_term(number_1)%j), &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(1),  &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_2_term(number_2)%j), &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2),  &
                                'coeff:',all_expansions%               &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%coeffs(icsf)
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                              'jj2 coupling' .and.                     &
                              all_expansions%coupling_expansions(      &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 2) then
                                number_1 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(2),  &
                                2,1000)
                                number_2 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2),  &
                                2,1000)
                                number_3 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(2),  &
                                1,1000)
                                number_4 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2),  &
                                1,1000)
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),2,1000),jj_1_term,num_1)
                               else
                                  jj_1_term(number_1)%j = -1
                                  jj_1_term(number_1)%subshellJ = 0
                                  jj_1_term(number_1)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),2,1000),jj_2_term,num_2)
!
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),1,1000),jj_3_term,num_3)
                               else
                                  jj_3_term(number_3)%j = -1
                                  jj_3_term(number_3)%subshellJ = 0
                                  jj_3_term(number_3)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),1,1000),jj_4_term,num_4)
                                write(iwrite_expansions,               &
                                '(1x,i5,3x,i2,a2,"(",i1,")",           &
                                3(i2,a2,"(",i1,")",6x),a6,f10.7)')icsf,&
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_1_term(number_1)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),2,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_2_term(number_2)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),2,1000),                        &
!
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_3_term(number_3)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),1,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_4_term(number_4)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),1,1000),                        &
!
                                'coeff:',all_expansions%               &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%coeffs(icsf)

                           else if(all_expansions%coupling_expansions( &
                              icoupling)%expansions(istate)%csfs(icsf)%&
                                                       nosubc == 2) then
                              write(iwrite_expansions,                 &
                                '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,  &
                                f10.7)') icsf,(all_expansions%         &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%in,CVAL(1,all_expansions%  &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%il),all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%iN_big,j=1,2),'coeff:',    &
                                all_expansions%coupling_expansions(    &
                                icoupling)%expansions(istate)%         &
                                coeffs(icsf)
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                              'cLSJ3 coupling' .and.                   &
                              all_expansions%coupling_expansions(      &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                number_1 = all_expansions%             &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(3)
                                number_2 = all_expansions%             &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(3)
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il -1, all_expansions%   &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),jj_1_term,num_1)
                               else
                                  jj_1_term(number_1)%j = -1
                                  jj_1_term(number_1)%subshellJ = 0
                                  jj_1_term(number_1)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il+1,all_expansions%       &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(1),  &
                                jj_2_term,num_2)
                                write(iwrite_expansions,               &
                                '(1x,i5,3x,2(i2,a2,"(",i1,")"),        &
                                2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_1_term(number_1)%j), &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(1),  &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_2_term(number_2)%j), &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(1),  &
                                (all_expansions%                       &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%in,                        &
                                CVAL(1,all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%il),all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%iN_big, j=2,3),            &
                                'coeff:',all_expansions%               &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%coeffs(icsf)

!
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                              'jj3 coupling' .and.                     &
                              all_expansions%coupling_expansions(      &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                number_1 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(2),  &
                                2,1000)
                                number_2 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2),  &
                                2,1000)
                                number_3 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM1(2),  &
                                1,1000)
                                number_4 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(2),  &
                                1,1000)
                                number_5 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(3),  &
                                2,1000)
                                number_6 =                             &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iM2(3),  &
                                1,1000)
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(1)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),2,1000),jj_1_term,num_1)
                               else
                                  jj_1_term(number_1)%j = -1
                                  jj_1_term(number_1)%subshellJ = 0
                                  jj_1_term(number_1)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),2,1000),jj_2_term,num_2)
!
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(2)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(1),1,1000),jj_3_term,num_3)
                               else
                                  jj_3_term(number_3)%j = -1
                                  jj_3_term(number_3)%subshellJ = 0
                                  jj_3_term(number_3)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),1,1000),jj_4_term,num_4)
!
                               if(all_expansions%                      &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(3)%il /=0 ) then
                                  call gettermjj (2*all_expansions%    &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  subc_cfg(3)%il -1,                   &
                                  JTHN(all_expansions%                 &
                                  coupling_expansions(icoupling)%      &
                                  expansions(istate)%csfs(icsf)%       &
                                  iM1(3),2,1000),jj_5_term,num_5)
                               else
                                  jj_5_term(number_5)%j = -1
                                  jj_5_term(number_5)%subshellJ = 0
                                  jj_5_term(number_5)%nu  = 0
                               end if
                               call gettermjj(2*all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(3)%il+1,                      &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(3),1,1000),jj_6_term,num_6)
                                write(iwrite_expansions,               &
                                '(1x,i5,3x,i2,a2,"(",i1,")",           &
                                5(i2,a2,"(",i1,")",6x),a6,f10.7)')icsf,&
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_1_term(number_1)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),2,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(1)%il,jj_2_term(number_2)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),2,1000),                        &
!
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_3_term(number_3)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(1),1,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(2)%il,jj_4_term(number_4)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM2(1),1,1000),                        &
!
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(3)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(3)%il,jj_5_term(number_5)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(3),2,1000),                        &
                                all_expansions%                        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(3)%in,                        &
                                LSJVAL(all_expansions%                 &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(3)%il,jj_6_term(number_6)%j), &
                                JTHN(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                iM1(3),1,1000),                        &
!
                                'coeff:',all_expansions%               &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%coeffs(icsf)

                           else if(all_expansions%coupling_expansions( &
                              icoupling)%expansions(istate)%csfs(icsf)%&
                              nosubc == 3  .and.  I_Numb == 3) then
                              write(iwrite_expansions,                 &
                                '(1x,i5,3x,3(i2,a1,"(",i1,")"),6x,a6,  &
                                f10.7)') icsf,(all_expansions%         &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%in,CVAL(1,all_expansions%  &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%il),all_expansions%        &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%         &
                                subc_cfg(j)%iN_big,j=1,3),'coeff:',    &
                                all_expansions%coupling_expansions(    &
                                icoupling)%expansions(istate)%         &
                                coeffs(icsf)
                           else
                              write(*,*) 'To many coupled shells in  ',&
                                 'Coupling_main.f90 printas'
                              write(*,*) 'nosubc, I_Numb',             &
                              all_expansions%coupling_expansions(      &
                              icoupling)%expansions(istate)%csfs(icsf)%&
                              nosubc, I_Numb
                              stop
                           end if
                        else
                           write(iwrite_expansions,*)                  &
                           '     all_expansions%coupling_expansions(', &
                         ')%expansions()%csfs()%subc_cfg NOT associated'
                        end if
!......................................................................
                        if(associated(all_expansions%                  &
                          coupling_expansions(icoupling)%expansions(   &
                          istate)%csfs(icsf)%subc) .and. associated(   &
                          all_expansions%coupling_expansions(          &
                          icoupling)%expansions(istate)%csfs(icsf)%iM1)&
                          .and.associated(all_expansions%              &
                          coupling_expansions(icoupling)%expansions(   &
                                           istate)%csfs(icsf)%iM2)) then
                        if(associated(all_expansions%                  &
                          coupling_expansions(icoupling)%expansions(   &
                                                   istate)%coeffs)) then
                           if(coupling_descriptions(all_expansions%    &
                             coupling_expansions(icoupling)%icoupling)%&
                                        long_name == 'LS coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 1) then
                                 write(iwrite_expansions,              &
                                 '(13x,i2,a1,i1)')                     &
                                 all_expansions%coupling_expansions(   &
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iS+1,CVAL(2,         &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    1)%inr
                              else if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc  &
                                                              == 2) then
                                 write(iwrite_expansions,              &
                                    '(12x,7(1x,i2,a1,i1))')            &
                                    (all_expansions%coupling_expansions&
                                    (icoupling)%expansions(istate)%    &
                                    csfs(icsf)%subc(j)%iS+1,CVAL(2,    &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(j)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,2),all_expansions%      &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM2(2)+1,CVAL(2,all_expansions%    &
                                    coupling_expansions(icoupling)%    &
                                   expansions(istate)%csfs(icsf)%iM1(2))
                              else if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)% nosubc &
                                 == 3  .and.  I_Numb == 3) then
                                 write(iwrite_expansions,              &
                                    '(12x,3(1x,i2,a1,i1),2(1x,i2,a1))')&
                                    (all_expansions%coupling_expansions&
                                    (icoupling)%expansions(istate)%    &
                                    csfs(icsf)%subc(j)%iS+1,CVAL(2,    &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(j)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,3),(all_expansions%     &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM2(j)+1,CVAL(2,all_expansions%    &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM1(j)),j=2,3)
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                    'jj1 coupling') then
                              if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc  &
                                 == 1) then
                                 write(iwrite_expansions,              &
                                '(9x,2(a4,1x,i1,1x),"[",a4,"]")')      &
                                JVAL(jj_1_term(number_1)%subshellJ),   &
                                jj_1_term(number_1)%nu,                &
                                JVAL(jj_2_term(number_2)%subshellJ),   &
                                jj_2_term(number_2)%nu,                &
                                JVAL(states%states(istate)%J)
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                     'JJ coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 1) then
                                 write(iwrite_expansions,              &
                                    '(13x,i2,a1,i1)')                  &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iS+1,CVAL(2,         &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    1)%inr
                              else if(all_expansions%                  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                                       nosubc == 2) then
                                 write(iwrite_expansions,              &
                                    '(12x,2(1x,i2,a1,i1),3x,"(",a4,",",&
                                    a4,")",a4)')(all_expansions%       &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iS+1,CVAL(2,all_expansions%     &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,2),JVAL(all_expansions% &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    1)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    2)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2))
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                     'LK coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 1) then
                                 write(iwrite_expansions,              &
                                    '(13x,i2,a1,i1)')                  &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iS+1,CVAL(2,         &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    1)%inr
                              else if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc  &
                                 == 2) then
                                 write(iwrite_expansions,              &
                                    '(12x,2(1x,i2,a1,i1),3x,a1," [",a4,&
                                    "]",a4)') (all_expansions%         &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iS+1,CVAL(2,all_expansions%     &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,2),CVAL(2,              &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%iM1(2)),JVAL(all_expansions% &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2)),JVAL(states%states(istate)%J)
                              end if

                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                     'JK coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 1) then
                                 write(iwrite_expansions,              &
                                    '(13x,i2,a1,i1)')                  &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iS+1,CVAL(2,         &
                                    all_expansions%coupling_expansions(&
                                    icoupling)%expansions(istate)%csfs(&
                                    icsf)%subc(1)%iL),all_expansions%  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    1)%inr
                              else if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc  &
                                                              == 2) then
                                 write(iwrite_expansions,              &
                                    '(12x,2(1x,i2,a1,i1),3x,a4," [",a4,&
                                    "]",a4)') (all_expansions%         &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(j)%iS+1,CVAL(2,all_expansions%&
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,2),JVAL(all_expansions% &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    1)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    2)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2))
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                  'LScjj coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 2) then
                                 write(iwrite_expansions,              &
                                 '(13x,i2,a1,i1,                       &
                                 2(a4,i2,1x)," [",a4,","a4," ]")')     &
                                 all_expansions%                       &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%        &
                                 subc(1)%iS+1,CVAL(2,all_expansions%   &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%subc(1)%&
                                 iL),all_expansions%                   &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%subc(1)%&
                                 inr,                                  &
                                 JVAL(jj_1_term(number_1)%subshellJ),  &
                                 jj_1_term(number_1)%nu,               &
                                 JVAL(jj_2_term(number_2)%subshellJ),  &
                                 jj_2_term(number_2)%nu,               &
                                 JVAL(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iJ(1)), &
                                 JVAL(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iJ(2))
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                    'jj2 coupling') then
                              if(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc  &
                                 == 2) then
                                 write(iwrite_expansions,              &
                                '(9x,2(a4,1x,i1,1x),"[",a4,"]",        &
                                2(a4,1x,i1,1x,"[",a4,"]"))')           &
                                JVAL(jj_1_term(number_1)%subshellJ),   &
                                jj_1_term(number_1)%nu,                &
                                JVAL(jj_2_term(number_2)%subshellJ),   &
                                jj_2_term(number_2)%nu,                &
                                JVAL(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(1)),  &
                                JVAL(jj_3_term(number_3)%subshellJ),   &
                                jj_3_term(number_3)%nu,                &
                                JVAL(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(2)),  &
                                JVAL(jj_4_term(number_4)%subshellJ),   &
                                jj_4_term(number_4)%nu,                &
                                JVAL(states%states(istate)%J)
                              end if

                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                    'LS3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                    '(12x,i2,a1,i1," (",i2,a1,i1,1x,i2,&
                                    a1,i1,")",i2,a1,1x,i2,a1)')        &
                                    (all_expansions%                   &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(j)%iS+1,CVAL(2,all_expansions%&
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,3),                     &
                                    (all_expansions% &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    j)+1,CVAL(2,all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    j)), j=2,3)
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                   'LSJ3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                   '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x, &
                                   i2,a1,i1,")",i2,a1,1x," [",a4,","a4,&
                                   " ]")')                             &
                                    (all_expansions%                   &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(j)%iS+1,CVAL(2,all_expansions%&
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,3),                     &
                                    all_expansions%                    &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2)+1,CVAL(2,all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    2)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iJ(1)),JVAL(all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                   expansions(istate)%csfs(icsf)%iJ(2))
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                   'LK3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                   '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x, &
                                   i2,a1,i1,")",i2,a1,1x,a1," [",a4,   &
                                   " ]")')                             &
                                    (all_expansions%                   &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(j)%iS+1,CVAL(2,all_expansions%&
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,3),                     &
                                    all_expansions%                    &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2)+1,CVAL(2,all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    2)),CVAL(2,all_expansions%         &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM1(3)),JVAL(all_expansions%       &
                                    coupling_expansions(icoupling)%    &
                                  expansions(istate)%csfs(icsf)%iM2(3))
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                   'JK3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                   '(12x,i2,a1,i1,1x," (",i2,a1,i1,1x, &
                                   i2,a1,i1,")",i2,a1,1x,a4,";",1x,a4, &
                                   " [",a4," ]")')                     &
                                    (all_expansions%                   &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(j)%iS+1,CVAL(2,all_expansions%&
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%iL),all_expansions%             &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%subc(&
                                    j)%inr,j=1,3),                     &
                                    all_expansions%                    &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM2( &
                                    2)+1,CVAL(2,all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%iM1( &
                                    2)),JVAL(all_expansions%           &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iJ(1)),JVAL(all_expansions%        &
                                    coupling_expansions(icoupling)%    &
                                  expansions(istate)%csfs(icsf)%iJ(2)),&
                                  JVAL(states%states(istate)%J)
                              end if
                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                   'cLSJ3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                 '(10x,2(1x,a4,i2)," (",i2,a1,i1,1x,i2,&
                                 a1,i1,")",i2,a1,1x," [",a4,","a4,     &
                                 " ]")')                               &
                                 JVAL(jj_1_term(number_1)%subshellJ),  &
                                 jj_1_term(number_1)%nu,               &
                                 JVAL(jj_2_term(number_2)%subshellJ),  &
                                 jj_2_term(number_2)%nu,               &
                                 (all_expansions%                      &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%        &
                                 subc(j)%iS+1,CVAL(2,all_expansions%   &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%subc(j)%&
                                 iL),all_expansions%                   &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%subc(j)%&
                                 inr, j=2,3),all_expansions%           &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iM2(2)  &
                                 +1,CVAL(2,all_expansions%             &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iM1(2)),&
                                 JVAL(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iJ(1)), &
                                 JVAL(all_expansions%                  &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%iJ(2))
                              end if

                           else if(coupling_descriptions(              &
                              all_expansions%coupling_expansions(      &
                              icoupling)%icoupling)%long_name ==       &
                                                   'jj3 coupling') then
                              if(all_expansions%coupling_expansions(   &
                                 icoupling)%expansions(istate)%csfs(   &
                                                 icsf)%nosubc == 3) then
                                 write(iwrite_expansions,              &
                                '(9x,2(a4,1x,i1,1x),"[",a4,"]",        &
                                4(a4,1x,i1,1x,"[",a4,"]"))')           &
                                JVAL(jj_1_term(number_1)%subshellJ),   &
                                jj_1_term(number_1)%nu,                &
                                JVAL(jj_2_term(number_2)%subshellJ),   &
                                jj_2_term(number_2)%nu,                &
                                JVAL(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(1)),  &
                                JVAL(jj_3_term(number_3)%subshellJ),   &
                                jj_3_term(number_3)%nu,                &
                                JVAL(JTHN(all_expansions%              &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(2),   &
                                2,1000)),                              &
                                JVAL(jj_4_term(number_4)%subshellJ),   &
                                jj_4_term(number_4)%nu,                &
                                JVAL(JTHN(all_expansions%              &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(2),   &
                                1,1000)),                              &
                                JVAL(jj_5_term(number_5)%subshellJ),   &
                                jj_5_term(number_5)%nu,                &
                                JVAL(all_expansions%                   &
                                coupling_expansions(icoupling)%        &
                                expansions(istate)%csfs(icsf)%iJ(3)),  &
                                JVAL(jj_6_term(number_6)%subshellJ),   &
                                jj_6_term(number_6)%nu,                &
                                JVAL(states%states(istate)%J)
                              end if
                           else
                             write(iwrite_expansions,'(8x,a5)')'state'
                             write(iwrite_expansions,'(9x,a20,a3,a2,   &
                                a3)')'subc  nr   L    S   ',           &
                                coupling_descriptions(all_expansions%  &
                                coupling_expansions(icoupling)%        &
                                icoupling)%iM1_name,'  ',              &
                                coupling_descriptions(all_expansions%  &
                                coupling_expansions(icoupling)%        &
                                icoupling)%iM2_name
                              do isubc=1,all_expansions%               &
                                 coupling_expansions(icoupling)%       &
                                 expansions(istate)%csfs(icsf)%nosubc,1
                                 inr = all_expansions%                 &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(isubc)%inr
                                 Li = all_expansions%                  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(isubc)%iL
                                 Si = all_expansions%                  &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    subc(isubc)%iS
                                 S_i = all_expansions%                 &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM2(isubc)
                                 L_i = all_expansions%                 &
                                    coupling_expansions(icoupling)%    &
                                    expansions(istate)%csfs(icsf)%     &
                                    iM1(isubc)
                                 write(iwrite_expansions,'(9x,i2,3x,i2,&
                                    3x,i2,3x,i2,3x,i2,3x,i2)')isubc,   &
                                    inr,Li,Si,L_i,S_i
                              end do
                           end if
                        else
                           write(iwrite_expansions,*)                  &
                             '    all_expansions%coupling_expansions(',&
                             ')%expansions()%coeffs NOT associated'
                        end if
                        else
                           write(iwrite_expansions,*)                  &
                           '     all_expansions%coupling_expansions(', &
                           ')%expansions()%csfs()%subc, iM1, iM2 NOT', &
                           ' associated'
                        end if
                        write(iwrite_expansions,*)' '
                     end do
                  else
                     write(iwrite_expansions,*) '      ',              &
                     'all_expansions%coupling_expansions(',            &
                     ')%expansions()%csfs NOT associated'
                  end if
               end do
            else
               write(iwrite_expansions, *) &
               'all_expansions%coupling_expansions()%expansions',&
               ' NOT associated'
            end if
         end do
      else
         write(iwrite_expansions,*) 'expansions NOT allocated'
      end if
!      write(*,*) ' end subroutine printas'
!
      end subroutine printas
!
!***********************************************************************
!                                                                      *
      subroutine print_DataBasis(icase,thresh)
!                                                                      *
!     This subroutine prints the oneconfigurational                    *
!     expansions to the files *.lbl                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      implicit none
      real(kind=dp), intent(in) :: thresh
      integer, intent(in) :: icase
      character(len=30) :: string_tmp
      character(len=9) :: string_shell_1, string_shell_2, string_shell_3
      character(len=9) :: string_shell_4, string_shell_5, string_shell_6
      character(len=9) :: string_shell_11, string_shell_12
      character(len=9) :: string_shell_13, string_shell_14
      integer          :: num_1, num_2, number_1, number_2
      integer          :: num_3, num_4, number_3, number_4
      integer          :: num_5, num_6, number_5, number_6
      character(len=1) :: CVAL
      character(len=2) :: LSJVAL
      character(len=4) :: JVAL
!      character(len=4) :: string_case_1, string_case_2, string_case_3
      character(len=14) :: string_case_1, string_case_2, string_case_3
!      character(len=4) :: string_case_4
      character(len=14) :: string_case_4, string_case_5, string_case_6
      character(len=2) :: string_shellX1, string_shellX2
      real(kind=dp)    :: Coeff
      real(kind=dp), dimension (:), pointer :: Coeff_Ind
      integer, dimension (:), pointer       :: indeksas
      type(subshell_term), dimension(1:63)  :: jj_1_term, jj_2_term
      type(subshell_term), dimension(1:63)  :: jj_3_term, jj_4_term
      type(subshell_term), dimension(1:63)  :: jj_5_term, jj_6_term
      integer :: icoupling, istate, icsf, j, ii, string_length, i_R
      if(associated(all_expansions%coupling_expansions)) then
        do icoupling=1, all_expansions%nr_of_coupling_expansions , 1
          if(associated(all_expansions%coupling_expansions(            &
            icoupling)%expansions)) then
            if(all_expansions%coupling_expansions(icoupling)%          &
              nr_of_expansions.ne.states%nr_of_states) then
              write(*,*) 'stop at subroutine printas: ',               &
                         'all_expansions%coupling_expansions(',        &
                         ')%nr_of_expansions.ne.states%nr_of_states'
              stop
            end if
            do istate=1, states%nr_of_states
              if(associated(all_expansions%coupling_expansions(        &
                               icoupling)%expansions(istate)%csfs)) then
                allocate(                                              &
                Coeff_Ind(all_expansions%coupling_expansions(          &
                icoupling)%expansions(istate)%size+1))
                allocate(                                              &
                indeksas(all_expansions%coupling_expansions(           &
                icoupling)%expansions(istate)%size+1))
                Coeff = 0.0D00
                Coeff_Ind = 0.0D00
                indeksas = 0
                do icsf=1,all_expansions%coupling_expansions(          &
                                    icoupling)%expansions(istate)%size,1
                  if(associated(all_expansions%                        &
                        coupling_expansions(icoupling)%expansions(     &
                                    istate)%csfs(icsf)%subc_cfg)) then
                    Coeff = Coeff +                                    &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%coeffs(icsf)                    &
                    * all_expansions%coupling_expansions(icoupling)%   &
                    expansions(istate)%coeffs(icsf)
!
                    if(icsf == 1) then
                      Coeff_Ind(1) =                                   &
                      all_expansions%coupling_expansions(icoupling)%   &
                      expansions(istate)%coeffs(1)
                      indeksas(1) = 1
                    else
                      do j = 1,icsf
                        if(j == icsf) then
                          Coeff_Ind(icsf) = all_expansions%            &
                                       coupling_expansions(icoupling)% &
                                       expansions(istate)%coeffs(icsf)
                          indeksas(icsf) = icsf
                        else
                          if(dabs(all_expansions%                      &
                          coupling_expansions(icoupling)%              &
                          expansions(istate)%coeffs(icsf)) >=          &
                                               dabs(Coeff_Ind(j))) then
                            do ii = icsf,j,-1
                               Coeff_Ind(ii+1) = Coeff_Ind(ii)
                               indeksas(ii+1) = indeksas(ii)
                            end do
                            Coeff_Ind(j) = all_expansions%             &
                                       coupling_expansions(icoupling)% &
                                       expansions(istate)%coeffs(icsf)
                            indeksas(j) = icsf
                            exit
                          end if
                        end if
                      end do
                    end if
                  end if
                end do
              end if
              do j = 1, all_expansions%coupling_expansions(            &
                                     icoupling)%expansions(istate)%size
                icsf = indeksas(j)
                if(dabs(all_expansions%coupling_expansions(icoupling)% &
                    expansions(istate)%coeffs(icsf)                    &
                    * all_expansions%coupling_expansions(icoupling)%   &
                    expansions(istate)%coeffs(icsf)) < dabs(thresh))exit
                if(coupling_descriptions(all_expansions%               &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj1 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 1) then
                  number_1=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(1),2,1000)
                  number_2=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(1),1,1000)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term,num_1)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il,     &
                     jj_1_term(number_1)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term(number_1)%subshellJ,                    &
                     jj_1_term(number_1)%nu,string_shell_1)
!                     jj_1_term(number_1)%Nr,string_shell_1)
                  else
                     jj_1_term(number_1)%j = -1
                     jj_1_term(number_1)%subshellJ = 0
                     jj_1_term(number_1)%nu  = 0
                     string_shell_1 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM1(1),1,1000),        &
                  jj_2_term,num_2)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  jj_2_term(number_2)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM1(1),1,1000),        &
                  jj_2_term(number_2)%subshellJ,                       &
                  jj_2_term(number_2)%nu,string_shell_2)
!                  jj_2_term(number_2)%Nr,string_shell_2)
!
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj2 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  number_1=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(2),2,1000)
                  number_2=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(2),2,1000)
                  number_3=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(2),1,1000)
                  number_4=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(2),1,1000)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term,num_1)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il,     &
                     jj_1_term(number_1)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term(number_1)%subshellJ,                    &
                     jj_1_term(number_1)%nu,string_shell_1)
!                     jj_1_term(number_1)%Nr,string_shell_1)
                  else
                     jj_1_term(number_1)%j = -1
                     jj_1_term(number_1)%subshellJ = 0
                     jj_1_term(number_1)%nu  = 0
                     string_shell_1 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),2,1000),        &
                  jj_2_term,num_2)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  jj_2_term(number_2)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),2,1000),        &
                  jj_2_term(number_2)%subshellJ,                       &
                  jj_2_term(number_2)%nu,string_shell_2)
!                  jj_2_term(number_2)%Nr,string_shell_2)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),1,1000),     &
                     jj_3_term,num_3)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il,     &
                     jj_3_term(number_3)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),1,1000),     &
                     jj_3_term(number_3)%subshellJ,                    &
                     jj_3_term(number_3)%nu,string_shell_3)
!                     jj_3_term(number_3)%Nr,string_shell_3)
                  else
                     jj_3_term(number_3)%j = -1
                     jj_3_term(number_3)%subshellJ = 0
                     jj_3_term(number_3)%nu  = 0
                     string_shell_3 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),1,1000),        &
                  jj_4_term,num_4)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il,        &
                  jj_4_term(number_4)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),1,1000),        &
                  jj_4_term(number_4)%subshellJ,                       &
                  jj_4_term(number_4)%nu,string_shell_4)
!                  jj_4_term(number_4)%Nr,string_shell_4)
!
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'cLSJ3 coupling' .and.                                 &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  number_1=                                            &
                         all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(3)
                  number_2=                                            &
                         all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(3)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il -1,  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%iM1(1),             &
                     jj_1_term,num_1)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il,     &
                     jj_1_term(number_1)%j,                            &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%iM1(1),             &
                     jj_1_term(number_1)%subshellJ,                    &
                     jj_1_term(number_1)%nu,string_shell_1)
!                     jj_1_term(number_1)%Nr,string_shell_1)
                  else
                     jj_1_term(number_1)%j = -1
                     jj_1_term(number_1)%subshellJ = 0
                     jj_1_term(number_1)%nu  = 0
                     string_shell_1 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il +1,     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(1),                &
                  jj_2_term,num_2)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  jj_2_term(number_2)%j,                               &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(1),                &
                  jj_2_term(number_2)%subshellJ,                       &
                  jj_2_term(number_2)%nu,string_shell_2)
!                  jj_2_term(number_2)%Nr,string_shell_2)
                  call getchLS (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%iN_big,    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(2)%iL,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(2)%iS,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(2)%inr,           &
                  string_shell_3)
                  call getchLS (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%il,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%iN_big,    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(3)%iL,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(3)%iS,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(3)%inr,           &
                  string_shell_4)
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LScjj coupling' .and.                                 &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  number_1=                                            &
                         all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(2)
                  number_2=                                            &
                         all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(2)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il -1,  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%iM1(1),             &
                     jj_1_term,num_1)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il,     &
                     jj_1_term(number_1)%j,                            &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%iM1(1),             &
                     jj_1_term(number_1)%subshellJ,                    &
                     jj_1_term(number_1)%nu,string_shell_2)
!                     jj_1_term(number_1)%Nr,string_shell_2)
                  else
                     jj_1_term(number_1)%j = -1
                     jj_1_term(number_1)%subshellJ = 0
                     jj_1_term(number_1)%nu  = 0
                     string_shell_2 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il +1,     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(1),                &
                  jj_2_term,num_2)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il,        &
                  jj_2_term(number_2)%j,                               &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(1),                &
                  jj_2_term(number_2)%subshellJ,                       &
                  jj_2_term(number_2)%nu,string_shell_3)
!                  jj_2_term(number_2)%Nr,string_shell_3)
                  call getchLS (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%iN_big,    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%iL,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%iS,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%inr,           &
                  string_shell_1)
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj3 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  number_1=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(2),2,1000)
                  number_2=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(2),2,1000)
                  number_3=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM1(2),1,1000)
                  number_4=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(2),1,1000)
                  number_5=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(3),2,1000)
                  number_6=                                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%iM2(3),1,1000)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term,num_1)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%il,     &
                     jj_1_term(number_1)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),2,1000),     &
                     jj_1_term(number_1)%subshellJ,                    &
                     jj_1_term(number_1)%nu,string_shell_1)
!                     jj_1_term(number_1)%Nr,string_shell_1)
                  else
                     jj_1_term(number_1)%j = -1
                     jj_1_term(number_1)%subshellJ = 0
                     jj_1_term(number_1)%nu  = 0
                     string_shell_1 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),2,1000),        &
                  jj_2_term,num_2)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  jj_2_term(number_2)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),2,1000),        &
                  jj_2_term(number_2)%subshellJ,                       &
                  jj_2_term(number_2)%nu,string_shell_2)
!                  jj_2_term(number_2)%Nr,string_shell_2)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),1,1000),     &
                     jj_3_term,num_3)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il,     &
                     jj_3_term(number_3)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(1),1,1000),     &
                     jj_3_term(number_3)%subshellJ,                    &
                     jj_3_term(number_3)%nu,string_shell_3)
!                     jj_3_term(number_3)%Nr,string_shell_3)
                  else
                     jj_3_term(number_3)%j = -1
                     jj_3_term(number_3)%subshellJ = 0
                     jj_3_term(number_3)%nu  = 0
                     string_shell_3 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),1,1000),        &
                  jj_4_term,num_4)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%il,        &
                  jj_4_term(number_4)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM2(1),1,1000),        &
                  jj_4_term(number_4)%subshellJ,                       &
                  jj_4_term(number_4)%nu,string_shell_4)
!                  jj_4_term(number_4)%Nr,string_shell_4)
                  if(all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%il /=0 )&
                                                                    then
                     call gettermjj                                    &
                     (2*all_expansions%coupling_expansions(icoupling)% &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%il -1,  &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(3),2,1000),     &
                     jj_5_term,num_5)
                     call getchjj (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%il,     &
                     jj_5_term(number_5)%j,                            &
                    JTHN(all_expansions%coupling_expansions(icoupling)%&
                     expansions(istate)%csfs(icsf)%iM1(3),2,1000),     &
                     jj_5_term(number_5)%subshellJ,                    &
                     jj_5_term(number_5)%nu,string_shell_5)
!                     jj_5_term(number_5)%Nr,string_shell_5)
                  else
                     jj_5_term(number_5)%j = -1
                     jj_5_term(number_5)%subshellJ = 0
                     jj_5_term(number_5)%nu  = 0
                     string_shell_5 = ""
                  end if
                  call gettermjj                                       &
                  (2*all_expansions%coupling_expansions(icoupling)%    &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%il +1,     &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM1(3),1,1000),        &
                  jj_6_term,num_6)
                  call getchjj (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(3)%il,        &
                  jj_6_term(number_6)%j,                               &
                  JTHN(all_expansions%coupling_expansions(icoupling)%  &
                  expansions(istate)%csfs(icsf)%iM1(3),1,1000),        &
                  jj_6_term(number_6)%subshellJ,                       &
                  jj_6_term(number_6)%nu,string_shell_6)
!                  jj_6_term(number_6)%Nr,string_shell_6)
!
                else
                  call getchLS (                                       &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%il,        &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%iN_big,    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%iL,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%iS,            &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc(1)%inr,           &
                  string_shell_1)
!
                  if(all_expansions%coupling_expansions(icoupling)%    &
                         expansions(istate)%csfs(icsf)%nosubc >= 2) then
                     call getchLS (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%il,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%iN_big, &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(2)%iL,         &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(2)%iS,         &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(2)%inr,        &
                     string_shell_2)
                  end if
                  if(all_expansions%coupling_expansions(icoupling)%    &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                     call getchLS (                                    &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%il,     &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%iN_big, &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(3)%iL,         &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(3)%iS,         &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc(3)%inr,        &
                     string_shell_3)
                  end if
                end if
! .....................................................................
                if(coupling_descriptions(all_expansions%               &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj1 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 1) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_jj123,'(A53)')                 &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_jj123,               &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  if(                                                  &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)   &
                                                                   then
                    if(len_trim(trim(string_shell_1)) == 0) then
                       string_case_1 = ''
                    else if(string_shell_1(4:4) == '_') then
                       string_case_1 = ''
                    else
                       string_case_1 = '('//                           &
                     trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))&
                       //') '
                    end if
                  else if(                                             &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)   &
                                                                   then
                    if(len_trim(trim(string_shell_1)) == 0) then
                       string_case_1 = ''
                    else if(string_shell_1(5:5) == '_') then
                       string_case_1 = ''
                    else
                       string_case_1 = '('//                           &
                     trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))&
                       //') '
                    end if
                  end if
                  if(                                                  &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)   &
                                                                   then
                    if(len_trim(trim(string_shell_2)) == 0) then
                       string_case_2 = ''
                    else if(string_shell_2(4:4) == '_') then
                       string_case_2 = ''
                    else if(string_shell_2(5:5) == '_') then
                       string_case_2 = ''
                    else
                       string_case_2 = '('//                           &
                     trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))&
                       //')'
                    end if
                  else if(                                             &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)   &
                                                                   then
                    if(len_trim(trim(string_shell_2)) == 0) then
                       string_case_2 = ''
                    else if(string_shell_2(5:5) == '_') then
                       string_case_2 = ''
                    else if(string_shell_2(6:6) == '_') then
                       string_case_2 = ''
                    else
                       string_case_2 = '('//                           &
                     trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))&
                       //')'
                    end if
                  end if
!
                  write(iwrite_DataBase_jj123,                         &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_shell_1)//trim(string_case_1)//          &
                  trim(string_shell_2)//trim(string_case_2)//          &
                  ' <'//                                               &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LS coupling') then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LS,'(A53)')                    &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LS,                  &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  if(all_expansions%coupling_expansions(icoupling)%    &
                         expansions(istate)%csfs(icsf)%nosubc == 1) then
                     write(iwrite_DataBase_LS,                         &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//'<'//                       &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  else if(all_expansions%coupling_expansions(icoupling)%&
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                    call getchXLS (                                    &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM1(2),              &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM2(2),              &
                    string_shellX1)
                    if(len_trim(trim(string_shell_2)) == 2) then
                        string_tmp=trim(adjustl(string_shell_2(1:2)))  &
                        //'_ '//trim(string_shellX1)
                    else if(len_trim(trim(string_shell_2)) == 3        &
                               .and. string_shell_2(3:3) == "_") then
                        string_tmp=trim(adjustl(string_shell_2(1:3)))  &
                        //' '//trim(string_shellX1)
                    else if(len_trim(trim(string_shell_2)) == 3) then
                        string_tmp=trim(adjustl(string_shell_2(1:3)))  &
                        //'_ '//trim(string_shellX1)
                    else
                        string_tmp=trim(adjustl(string_shell_2))//" "//&
                        trim(string_shellX1)
                    end if
                    if(                                                &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)&
                                                                   then
                       i_R = 4
                    else if(                                           &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)&
                                                                   then
                       i_R = 5
                    end if
                    if(len_trim(trim(string_shell_1)) <= i_R) then
                       write(iwrite_DataBase_LS,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1)//                          &
                       trim(string_tmp)//'<'//                         &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                    else
                       write(iwrite_DataBase_LS,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1)//'.'//                     &
                       trim(string_tmp)//'<'//                         &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                     end if

                  else if(all_expansions%coupling_expansions(          &
                    icoupling)%expansions(istate)%csfs(icsf)%nosubc    &
                                                              == 3) then
                    call getchXLS (                                    &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM1(2),              &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM2(2),              &
                    string_shellX1)
                    call getchXLS (                                    &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM1(3),              &
                    all_expansions%coupling_expansions(icoupling)%     &
                    expansions(istate)%csfs(icsf)%iM2(3),              &
                    string_shellX2)
!
                    if(len_trim(trim(string_shell_2)) == 3) then
                       write(iwrite_DataBase_LS,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1)//'.'//                     &
                       trim(string_shell_2)//                          &
                       trim(string_shellX1)//'.'//                     &
                       trim(string_shell_3)//' '//                     &
                       string_shellX2//'<'//                           &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                    else
                       write(iwrite_DataBase_LS,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1)//'.'//                     &
                       trim(string_shell_2)//                          &
                       trim(string_shellX1)//'.'//                     &
                       trim(string_shell_3)//' '//                     &
                       string_shellX2//'<'//                           &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                     end if
                   end if
!
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'JJ coupling') then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_JJ,'(A53)')                    &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_JJ,                  &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  if(all_expansions%coupling_expansions(icoupling)%    &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                    if(len_trim(trim(string_shell_2)) == 3) then
                       string_tmp=trim(adjustl(string_shell_2(1:2)))   &
                       //'<'//trim(adjustl                             &
                       (JVAL(all_expansions%coupling_expansions(       &
                       icoupling)%expansions(istate)%csfs(icsf)%       &
                       iM1(2))))//'>'
                    else if(len_trim(trim(string_shell_2)) == 4) then
                       string_tmp=trim(adjustl(string_shell_2(1:3)))   &
                       //'<'//trim(adjustl                             &
                       (JVAL(all_expansions%coupling_expansions(       &
                       icoupling)%expansions(istate)%csfs(icsf)%       &
                       iM1(2))))//'>'
                    else
                       string_length =                                 &
                               len_trim(trim(adjustl(string_shell_2)))-1
                       string_tmp=trim(adjustl(string_shell_2          &
                       (1:string_length)))//'<'//trim(adjustl          &
                       (JVAL(all_expansions%coupling_expansions(       &
                       icoupling)%expansions(istate)%csfs(icsf)%       &
                       iM1(2))))//'>)'
                    end if
                    string_length =                                    &
                               len_trim(trim(adjustl(string_shell_1)))-1
                    if(len_trim(trim(string_shell_1)) <= 4) then
                       write(iwrite_DataBase_JJ,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1(1:string_length))//'<'//    &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM2(1))))//'>.'//            &
                       trim(string_tmp)//' ('//                        &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM2(1))))//','//             &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM1(2))))//')<'//            &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                    else
                       string_length = len_trim(trim(string_shell_1))-1
                       write(iwrite_DataBase_JJ,                       &
                       '(7X,F12.8,3X,F11.8,3X,A)')                     &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf) *               &
                       all_expansions%coupling_expansions(icoupling)%  &
                       expansions(istate)%coeffs(icsf),                &
                       trim(string_shell_1(1:string_length))//'<'//    &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM2(1))))//'>).'//           &
                       trim(string_tmp)//' ('//                        &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM2(1))))//','//             &
                       trim(adjustl(JVAL(all_expansions%               &
                       coupling_expansions(icoupling)%expansions(      &
                       istate)%csfs(icsf)%iM1(2))))//')<'//            &
                       trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                     end if
                   end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LK coupling' .and.                                    &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LK,'(A53)')                    &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LK,                  &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*101,"%"
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9 &
                     .and.                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) /= 3 ) then
                      string_tmp= string_shell_2(5:5)
                    else if(all_expansions%coupling_expansions(        &
                      icoupling)%expansions(istate)%csfs(icsf)%subc(2)%&
                       iS /= 0) then
                       string_tmp= adjustl(JVAL(2*all_expansions%      &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%subc(2)%iS+2))
                    else if(len_trim(trim(string_shell_1)) /= 3 ) then
                       string_tmp= string_shell_1(5:5)
                    else
                       string_tmp= '1'
                    end if
                    i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9 &
                     .and.                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) /= 4 ) then
                      string_tmp= string_shell_2(6:6)
                    else if(all_expansions%coupling_expansions(        &
                      icoupling)%expansions(istate)%csfs(icsf)%subc(2)%&
                       iS /= 0) then
                       string_tmp= adjustl(JVAL(2*all_expansions%      &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%subc(2)%iS+2))
                    else if(len_trim(trim(string_shell_1)) /= 3 ) then
                       string_tmp= string_shell_1(5:5)
                    else
                       string_tmp= '1'
                    end if
                    i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10 &
                     .and.                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) /= 3 ) then
                      string_tmp= string_shell_2(5:5)
                    else if(all_expansions%coupling_expansions(        &
                      icoupling)%expansions(istate)%csfs(icsf)%subc(2)%&
                       iS /= 0) then
                       string_tmp= adjustl(JVAL(2*all_expansions%      &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%subc(2)%iS+2))
                    else if(len_trim(trim(string_shell_1)) /= 4 ) then
                       string_tmp= string_shell_1(6:6)
                    else
                       string_tmp= '1'
                    end if
                    i_R = 5
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10 &
                     .and.                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) /= 4 ) then
                      string_tmp= string_shell_2(6:6)
                    else if(all_expansions%coupling_expansions(        &
                      icoupling)%expansions(istate)%csfs(icsf)%subc(2)%&
                       iS /= 0) then
                       string_tmp= adjustl(JVAL(2*all_expansions%      &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%subc(2)%iS+2))
                    else if(len_trim(trim(string_shell_1)) /= 4 ) then
                       string_tmp= string_shell_1(6:6)
                    else
                       string_tmp= '1'
                    end if
                    i_R = 5
                  end if
!
                  if(len_trim(trim(string_shell_1)) <= i_R) then
                     write(iwrite_DataBase_LK,                         &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//' '//                       &
                     trim(CVAL(2,all_expansions%coupling_expansions(   &
                     icoupling)%expansions(istate)%csfs(icsf)%iM1(2))) &
                     //'_'//trim(string_tmp)//'['//                    &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(                              &
                     icoupling)%expansions(istate)%csfs(icsf)%iM2(2))))&
                     //']<'//                                          &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  else
                     write(iwrite_DataBase_LK,                         &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//'.'//                       &
                     trim(string_shell_2)//' '//                       &
                     trim(CVAL(2,all_expansions%coupling_expansions(   &
                     icoupling)%expansions(istate)%csfs(icsf)%iM1(2))) &
                     //'_'//trim(string_tmp)//'['//                    &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(                              &
                     icoupling)%expansions(istate)%csfs(icsf)%iM2(2))))&
                     //']<'//                                          &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'JK coupling' .and.                                    &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_JK,'(A53)')                    &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_JK,                  &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)&
                                                                   then
                    if(len_trim(trim(string_shell_1)) == 3) then
                       string_tmp= string_shell_1(1:2)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>'
                    else if(len_trim(trim(string_shell_1)) == 4 ) then
                       string_tmp= string_shell_1(1:3)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>'
                    else
                       string_tmp= string_shell_1(1:6)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>)'
                    end if
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)&
                                                                   then
                    if(len_trim(trim(string_shell_1)) == 4) then
                       string_tmp= string_shell_1(1:3)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>'
                    else if(len_trim(trim(string_shell_1)) == 5 ) then
                       string_tmp= string_shell_1(1:4)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>'
                    else
                       string_tmp= string_shell_1(1:7)//'<'//          &
                       trim(adjustl(JVAL(all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                       csfs(icsf)%iM2(1))))//'>)'
                    end if
                  end if
!
!                  if(len_trim(trim(string_shell_2)) /= 3) then
!                     write(iwrite_DataBase_JK,                         &
!                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
!                     all_expansions%coupling_expansions(icoupling)%    &
!                     expansions(istate)%coeffs(icsf),                  &
!                     all_expansions%coupling_expansions(icoupling)%    &
!                     expansions(istate)%coeffs(icsf) *                 &
!                     all_expansions%coupling_expansions(icoupling)%    &
!                     expansions(istate)%coeffs(icsf),                  &
!                     trim(string_tmp)//'.'//                           &
!                     trim(string_shell_2)//'  '//                      &
!                     string_shell_2(5:5)//'['//                        &
!                     trim(adjustl(JVAL(all_expansions%                 &
!                     coupling_expansions(                              &
!                     icoupling)%expansions(istate)%csfs(icsf)%iM1(2))))&
!                     //']<'//                                          &
!                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
!                  else
                     write(iwrite_DataBase_JK,                         &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_tmp)//'.'//                           &
                     trim(string_shell_2)//' '//                       &
                     trim(adjustl(JVAL(2*all_expansions%               &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%subc(2)%iS+2)))//'['//                 &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(                              &
                     icoupling)%expansions(istate)%csfs(icsf)%iM1(2))))&
                     //']<'//                                          &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
!                  end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LScjj coupling' .and.                                 &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LScjj,'(A53)')                 &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LScjj,               &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(1),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(1),                &
                  string_shellX1)
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) == 0) then
                       string_case_2 = ''
                    else if(string_shell_2(4:4) /= '_' .and.           &
                               len_trim(trim(string_shell_3)) /= 0) then
                       string_case_2 = '.'
                    else
                       string_case_2 = ''
                    end if
                    if(string_shell_3(4:4) == '_') then
                       string_case_3 = '('
                    else if(string_shell_2(4:4) == '_' .and.           &
                                 len_trim(trim(string_shell_3)) == 0) then
                       string_case_3 = '('
                    else
                       string_case_3 = ' ('
                    end if
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                    if(len_trim(trim(string_shell_2)) == 0) then
                       string_case_2 = ''
                    else if(string_shell_2(5:5) /= '_' .and.           &
                               len_trim(trim(string_shell_3)) /= 0) then
                       string_case_2 = '.'
                    else
                       string_case_2 = ''
                    end if
                    if(string_shell_3(5:5) == '_') then
                       string_case_3 = '('
                    else if(string_shell_2(5:5) == '_' .and.           &
                                 len_trim(trim(string_shell_3)) == 0) then
                       string_case_3 = '('
                    else
                       string_case_3 = ' ('
                    end if
                  end if
                  write(iwrite_DataBase_LScjj,                         &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_shell_1)//'<'//                          &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(1))))//'> '//                          &
                  trim(string_shell_2)//trim(string_case_2)//          &
                  trim(string_shell_3)//trim(string_case_3)//          &
                  trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))// &
                  ','//                                                &
                  trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))// &
                  ')<'//                                               &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(2))))//'> ('//                         &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(1))))//','//                           &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(2))))//')<'//                          &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj2 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 2) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_jj123,'(A53)')                 &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_jj123,               &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  if(len_trim(trim(string_shell_1)) == 0) then
                     string_case_1 = ''
                  else if(string_shell_1(4:4) == '_') then
                     string_case_1 = ''
                  else if(string_shell_1(5:5) == '_') then
                     string_case_1 = ''
                  else if(string_shell_1(6:6) == '_') then
                     string_case_1 = ''
                  else
                     string_case_1 = '('//                             &
                     trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))&
                     //') '
                  end if
                  if(len_trim(trim(string_shell_2)) == 0) then
                     string_case_2 = ''
                  else if(string_shell_2(4:4) == '_') then
                     string_case_2 = ''
                  else if(string_shell_2(5:5) == '_') then
                     string_case_2 = ''
                  else if(string_shell_2(6:6) == '_') then
                     string_case_2 = ''
                  else
                     string_case_2 = '('//                             &
                     trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))&
                     //')'
                  end if
                  if(len_trim(trim(string_shell_1)) == 0) then
                  else if(len_trim(trim(string_shell_2)) == 0) then
                  else if(len_trim(trim(string_shell_3)) /= 0 .or.     &
                     len_trim(trim(string_shell_4)) /= 0) then
                     string_case_2 = trim(string_case_2)//'<'//        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//'>.'
                  end if
                  if(len_trim(trim(string_shell_3)) == 0) then
                     string_case_3 = ''
                  else if(string_shell_3(4:4) == '_') then
                     string_case_3 = ''
                  else if(string_shell_3(5:5) == '_') then
                     string_case_3 = ''
                  else if(string_shell_3(6:6) == '_') then
                     string_case_3 = ''
                  else
                     string_case_3 = '('//                             &
                     trim(adjustl(JVAL(jj_3_term(number_3)%subshellJ)))&
                     //')'
                  end if
                  if(len_trim(trim(string_shell_3)) /= 0 .and.         &
                     len_trim(trim(string_shell_4)) /= 0) then
                     string_case_3 = trim(string_case_3)//'<'//        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(2))))//'>.'
                  end if
                  if(len_trim(trim(string_shell_4)) == 0) then
                     string_case_4 = ''
                  else if(string_shell_4(4:4) == '_') then
                     string_case_4 = ''
                  else if(string_shell_4(5:5) == '_') then
                     string_case_4 = ''
                  else if(string_shell_4(6:6) == '_') then
                     string_case_4 = ''
                  else
                     string_case_4 = '('//                             &
                    trim(adjustl(JVAL(jj_4_term(number_4)%subshellJ))) &
                     //') '
                  end if
!
                  write(iwrite_DataBase_jj123,                         &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_shell_1)//trim(string_case_1)//          &
                  trim(string_shell_2)//trim(string_case_2)//          &
                  trim(string_shell_3)//trim(string_case_3)//          &
                  trim(string_shell_4)//trim(string_case_4)//          &
                  ' <'//                                               &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LS3 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then

                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LS3,'(A53)')                   &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LS3,                 &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
!
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(2),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(2),                &
                  string_shellX1)
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(3),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(3),                &
                  string_shellX2)
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 4
                  end if
!
                  if(len_trim(trim(string_shell_2)) == i_R) then
!                  if(len_trim(trim(string_shell_2)) == 3) then
                     write(iwrite_DataBase_LS3,                        &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//                            &
                     trim(string_shell_3)//                            &
                     '('//trim(string_shellX1)//') '//                 &
                     string_shellX2//'<'//                             &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  else
                     write(iwrite_DataBase_LS3,                        &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//'.'//                       &
                     trim(string_shell_3)//                            &
                     '('//trim(string_shellX1)//') '//                 &
                     string_shellX2//'<'//                             &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LSJ3 coupling' .and.                                  &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LSJ3,'(A53)')                  &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LSJ3,                &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
!
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(2),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(2),                &
                  string_shellX1)
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(3),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(3),                &
                  string_shellX2)
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 4
                  end if
                  if(len_trim(trim(string_shell_2)) == i_R) then
!                  if(len_trim(trim(string_shell_2)) == 3) then
                     write(iwrite_DataBase_LSJ3,                       &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//                            &
                     trim(string_shell_3)//                            &
                     '('//trim(string_shellX1)//') ('//                &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//','//                        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(2))))//')<'//                       &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  else
                     write(iwrite_DataBase_LSJ3,                       &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//'.'//                       &
                     trim(string_shell_3)//                            &
                     '('//trim(string_shellX1)//') ('//                &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//','//                        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(2))))//')<'//                       &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'LK3 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_LK3,'(A53)')                   &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_LK3,                 &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
!
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(2),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(2),                &
                  string_shellX1)
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 4
                  end if
                  if(len_trim(trim(string_shell_2)) == i_R) then
!                  if(len_trim(trim(string_shell_2)) == 3) then
                     write(iwrite_DataBase_LK3,                        &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//                            &
                     trim(string_shell_3)//'('//                       &
                     trim(string_shellX1)//') '//                      &
                     trim(CVAL(2,all_expansions%coupling_expansions(   &
                     icoupling)%expansions(istate)%csfs(icsf)%iM1(3))) &
                     //'_'//string_shellX1(1:1)//'['//                 &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(                              &
                     icoupling)%expansions(istate)%csfs(icsf)%iM2(3))))&
                     //']<'//                                          &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  else
                     write(iwrite_DataBase_LK3,                        &
                     '(7X,F12.8,3X,F11.8,3X,A)')                       &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf) *                 &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%coeffs(icsf),                  &
                     trim(string_shell_1)//' '//                       &
                     trim(string_shell_2)//'.'//                       &
                     trim(string_shell_3)//'('//                       &
                     trim(string_shellX1)//') '//                      &
                     trim(CVAL(2,all_expansions%coupling_expansions(   &
                     icoupling)%expansions(istate)%csfs(icsf)%iM1(3))) &
                     //'_'//string_shellX1(1:1)//'['//                 &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(                              &
                     icoupling)%expansions(istate)%csfs(icsf)%iM2(3))))&
                     //']<'//                                          &
                     trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                  end if
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'JK3 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_JK3,'(A53)')                   &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_JK3,                 &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
!
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(2),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(2),                &
                  string_shellX1)
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)&
                                                                   then
                       i_R = 4
                  end if
                  if(len_trim(trim(string_shell_1)) == i_R) then
!                  if(len_trim(trim(string_shell_1)) == 3) then
                     string_tmp= string_shell_1(1:2)//'<'//           &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//'>'
                  else if(string_shell_1(4:4) == '_') then
                     string_tmp= string_shell_1(1:4)//'<'//           &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//'>'
                  else
                     string_tmp= string_shell_1(1:6)//'<'//            &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//'>)'
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 4
                  end if
                  if(len_trim(trim(string_shell_2)) == i_R) then
!                  if(len_trim(trim(string_shell_2)) == 3) then
                     string_case_1 = ''
                  else
                     string_case_1 = '.'
                  end if
                  write(iwrite_DataBase_JK3,                           &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_tmp)//                                   &
                  trim(string_shell_2)//trim(string_case_1)//          &
                  trim(string_shell_3)//                               &
                  '('//trim(string_shellX1)//') '//                    &
                  string_shellX1(1:1)//'['//                           &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(                                 &
                  icoupling)%expansions(istate)%csfs(icsf)%iJ(2))))    &
                  //']<'//                                             &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'cLSJ3 coupling' .and.                                 &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_cLSJ3,'(A53)')                 &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_cLSJ3,               &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
                  call getchXLS (                                      &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM1(2),                &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%csfs(icsf)%iM2(2),                &
                  string_shellX1)
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_1)) == 0) then
                     string_case_1 = ''
!                  else if(string_shell_1(4:4) /= '_' .and.             &
                  else if(string_shell_1(i_R:i_R) /= '_' .and.         &
                               len_trim(trim(string_shell_2)) /= 0) then
                     string_case_1 = '.'
                  else
                     string_case_1 = ''
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(string_shell_2(i_R:i_R) == '_') then
                     string_case_2 = '('
                  else
                     string_case_2 = ' ('
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in <= 9)&
                                                                   then
                       i_R = 3
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in >=10)&
                                                                   then
                       i_R = 4
                  end if

                  if(len_trim(trim(string_shell_3)) == i_R) then
!                  if(len_trim(trim(string_shell_3)) == 3) then
                     string_case_3 = ''
                  else
                     string_case_3 = '.'
                  end if
                  write(iwrite_DataBase_cLSJ3,                         &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_shell_1)//trim(string_case_1)//          &
                  trim(string_shell_2)//trim(string_case_2)//          &
                  trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))// &
                  ','//                                                &
                  trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))// &
                  ')<'//                                               &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(1))))//'> '//                          &
                  trim(string_shell_3)//trim(string_case_3)//          &
                  trim(string_shell_4)//'('//trim(string_shellX1)//')' &
                  //'<'//                                              &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(2))))//'> ('//                         &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(1))))//','//                           &
                  trim(adjustl(JVAL(all_expansions%                    &
                  coupling_expansions(icoupling)%expansions(istate)%   &
                  csfs(icsf)%iJ(2))))//')<'//                          &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
!
                else if(coupling_descriptions(all_expansions%          &
                coupling_expansions(icoupling)%icoupling)%long_name == &
                'jj3 coupling' .and.                                   &
                all_expansions%coupling_expansions(icoupling)%         &
                         expansions(istate)%csfs(icsf)%nosubc == 3) then
                  if(istate == 1 .and. j ==1 .and. icase == 1)         &
                  write(iwrite_DataBase_jj123,'(A53)')                 &
                 " Pos   J   Parity      Energy Total      Comp. of ASF"
                  if(j ==1) write(iwrite_DataBase_jj123,               &
                  '(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)') istate,    &
                  trim(JVAL(states%states(istate)%J)),                 &
                  " ", states%states(istate)%energy,Coeff*100,"%"
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(1)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_1)) == 0) then
                     string_case_1 = ''
                  else if(string_shell_1(i_R:i_R) == '_') then
!                  else if(string_shell_1(4:4) == '_') then
                     string_case_1 = ''
                  else if(string_shell_1(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_1(5:5) == '_') then
                     string_case_1 = ''
                  else
                     string_case_1 = '('//                             &
                     trim(adjustl(JVAL(jj_1_term(number_1)%subshellJ)))&
                     //') '
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(2)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_2)) == 0) then
                     string_case_2 = ''
                  else if(string_shell_2(i_R:i_R) == '_') then
!                  else if(string_shell_2(4:4) == '_') then
                     string_case_2 = ''
                  else if(string_shell_2(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_2(5:5) == '_') then
                     string_case_2 = ''
                  else
                     string_case_2 = '('//                             &
                     trim(adjustl(JVAL(jj_2_term(number_2)%subshellJ)))&
                     //')'
                  end if
                  if(len_trim(trim(string_shell_1)) == 0) then
                  else if(len_trim(trim(string_shell_2)) == 0) then
                  else if(len_trim(trim(string_shell_3)) /= 0 .or.     &
                          len_trim(trim(string_shell_4)) /= 0 .or.     &
                          len_trim(trim(string_shell_5)) /= 0 .or.     &
                          len_trim(trim(string_shell_6)) /= 0) then
                     string_case_2 = trim(string_case_2)//'<'//        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(1))))//'>.'
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(3)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_3)) == 0) then
                     string_case_3 = ''
                  else if(string_shell_3(i_R:i_R) == '_') then
!                  else if(string_shell_3(4:4) == '_') then
                     string_case_3 = ''
                  else if(string_shell_3(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_3(5:5) == '_') then
                     string_case_3 = ''
                  else
                     string_case_3 = '('//                             &
                     trim(adjustl(JVAL(jj_3_term(number_3)%subshellJ)))&
                     //')'
                  end if
                  if(len_trim(trim(string_shell_3)) /= 0) then
                     if(len_trim(trim(string_shell_4)) /= 0 .or.       &
                        len_trim(trim(string_shell_5)) /= 0 .or.       &
                        len_trim(trim(string_shell_6)) /= 0) then
                        string_case_3 = trim(string_case_3)//'<'//     &
                        trim(adjustl(JVAL(JTHN(all_expansions%         &
                        coupling_expansions(icoupling)%expansions(     &
                        istate)%csfs(icsf)%iJ(2),2,1000))))//'>.'
                     end if
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(4)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(4)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_4)) == 0) then
                     string_case_4 = ''
                  else if(string_shell_4(i_R:i_R) == '_') then
!                  else if(string_shell_4(4:4) == '_') then
                     string_case_4 = ''
                  else if(string_shell_4(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_4(5:5) == '_') then
                     string_case_4 = ''
                  else
                     string_case_4 = '('//                             &
                    trim(adjustl(JVAL(jj_4_term(number_4)%subshellJ))) &
                     //') '
                  end if
                  if(len_trim(trim(string_shell_4)) /= 0) then
                     if(len_trim(trim(string_shell_5)) /= 0 .or.       &
                        len_trim(trim(string_shell_6)) /= 0) then
                        string_case_4 = trim(string_case_4)//'<'//     &
                        trim(adjustl(JVAL(JTHN(all_expansions%         &
                        coupling_expansions(icoupling)%expansions(     &
                        istate)%csfs(icsf)%iJ(2),1,1000))))//'>.'
                     end if
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(5)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(5)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_5)) == 0) then
                     string_case_5 = ''
                  else if(string_shell_5(i_R:i_R) == '_') then
!                  else if(string_shell_5(4:4) == '_') then
                     string_case_5 = ''
                  else if(string_shell_5(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_5(5:5) == '_') then
                     string_case_5 = ''
                  else
                     string_case_5 = '('//                             &
                    trim(adjustl(JVAL(jj_5_term(number_5)%subshellJ))) &
                     //') '
                  end if
                  if(len_trim(trim(string_shell_5)) /= 0 .and.         &
                     len_trim(trim(string_shell_6)) /= 0) then
                     string_case_5 = trim(string_case_5)//'<'//        &
                     trim(adjustl(JVAL(all_expansions%                 &
                     coupling_expansions(icoupling)%expansions(istate)%&
                     csfs(icsf)%iJ(3))))//'>.'
                  end if
!
                  if(                                                  &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(6)%in <= 9)&
                                                                   then
                       i_R = 4
                  else if(                                             &
                     all_expansions%coupling_expansions(icoupling)%    &
                     expansions(istate)%csfs(icsf)%subc_cfg(6)%in >=10)&
                                                                   then
                       i_R = 5
                  end if
                  if(len_trim(trim(string_shell_6)) == 0) then
                     string_case_6 = ''
                  else if(string_shell_6(i_R:i_R) == '_') then
!                  else if(string_shell_6(4:4) == '_') then
                     string_case_6 = ''
                  else if(string_shell_6(i_R+1:i_R+1) == '_') then
!                  else if(string_shell_6(5:5) == '_') then
                     string_case_6 = ''
                  else
                     string_case_6 = '('//                             &
                    trim(adjustl(JVAL(jj_6_term(number_6)%subshellJ))) &
                     //') '
                  end if
!
                  write(iwrite_DataBase_jj123,                         &
                  '(7X,F12.8,3X,F11.8,3X,A)')                          &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf) *                    &
                  all_expansions%coupling_expansions(icoupling)%       &
                  expansions(istate)%coeffs(icsf),                     &
                  trim(string_shell_1)//trim(string_case_1)//          &
                  trim(string_shell_2)//trim(string_case_2)//          &
                  trim(string_shell_3)//trim(string_case_3)//          &
                  trim(string_shell_4)//trim(string_case_4)//          &
                  trim(string_shell_5)//trim(string_case_5)//          &
                  trim(string_shell_6)//trim(string_case_6)//          &
                  ' <'//                                               &
                  trim(adjustl(JVAL(states%states(istate)%J)))//'>'
                end if
              end do
              deallocate(Coeff_Ind)
              deallocate(indeksas)
            end do
          end if
        end do
      else
        write(iwrite_expansions, *) &
        'all_expansions%coupling_expansions()%expansions',&
        ' NOT associated'
      end if
      end subroutine print_DataBasis
!
!---------------------------------------------------------------
!
      subroutine store_R_values(iwrite,note,couplings)
!--------------------------------------------------------------------
!     This subroutine stores the R_values of
!     the requested couplings to the unit "iwrite"
!
!--------------------------------------------------------------------
      implicit none
      integer, intent(in)::iwrite
      character(len=2), intent(in) :: note
      character(len=13), intent(in) :: couplings
      integer::nr_of_couplings
      integer :: j, jj, J_Max
      J_Max = len_trim(trim(couplings))
      jj = 0
      write(iwrite,*)                                                  &
      "-------------------------------------------------"
      write(iwrite,'(A)')                                              &
             ' List of  R  values in different coupling schemes:'
      write(iwrite,'(A)')                                              &
             ' (The smaller R the better is coupling)'
      do j =1,J_Max
        if(couplings(j:j).eq.'y') then
            jj = jj + 1
            write(iwrite,'(5x,a20,4x,f6.3)')                           &
            coupling_descriptions(j)%long_name,                        &
            all_evaluations%R(jj)%R_value
        end if
      end do
      end subroutine store_R_values
!
      subroutine store_P_values(iwrite,note,couplings)
!--------------------------------------------------------------------
!     This subroutine stores the P_values of
!     the requested couplings to the unit "iwrite"
!
!--------------------------------------------------------------------
      implicit none
      integer, intent(in) :: iwrite
      character(len=2), intent(in) :: note
      character(len=13), intent(in) :: couplings
      integer :: j, jj, J_Max
      J_Max = len_trim(trim(couplings))
      jj = 0
      write(iwrite,*)                                                  &
      "-------------------------------------------------"
      write(iwrite,'(A)')                                              &
             ' List of  P  values in different coupling schemes:'
      write(iwrite,'(A)')                                              &
             ' (The bigger P the better is coupling)'
      do j =1,J_Max
        if(couplings(j:j).eq.'y') then
            jj = jj + 1
            write(iwrite,'(5x,a20,4x,f6.3)')                           &
            coupling_descriptions(j)%long_name,                        &
            all_evaluations%R(jj)%P_value
        end if
      end do
      write(iwrite,*)                                                  &
      "-------------------------------------------------"
      end subroutine store_P_values
!
      end module Coupling_main
