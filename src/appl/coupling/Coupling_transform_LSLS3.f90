!
!***********************************************************************
!
      module Coupling_transform_LSLS3
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
      use Coupling_transform_cfg_LSLS3
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: main_LSLS3
      private :: cfg_count
!
      private :: form_cfg_structure
      private :: delete_cfg_structure
!-----------------------------------------------
!   D e f i n i t i o n s
!-----------------------------------------------

      type(cfg_structure), private::cfg_structure_LS
      integer, private :: nr_of_LS_coupling
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine main_LSLS3 (icoupling_nr_in_expansions,               &
                                         inr_of_LS_coupling,print_level)
!                                                                      *
!     This is managerial subroutine
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent (in):: icoupling_nr_in_expansions
      integer, intent (in):: inr_of_LS_coupling, print_level
      integer::istate
      integer::icfg, icsf_LS3
      integer::nr_of_csfs_LS3, cfg_nr_of_csfs_LS3
      real(kind=dp):: sum_of_state
!
!      write(*,*) '  subroutine main_lsls3'
      write(iwrite_log,*) ''
      write(iwrite_log,*) 'LS3 coupling:'
      nr_of_LS_coupling = inr_of_LS_coupling
      if(nr_of_LS_coupling.lt.1) then
         write(*,*)                                                    &
           'ERROR in subroutine main_LSLS3: nr_of_LS_coupling.lt.1; ', &
           'nr_of_LS_coupling=', nr_of_LS_coupling
         write(*,*)                                                    &
          '(not specified a serial number of expansions in LS coupling)'
         write(*,*)'program will be terminated'
         write(iwrite_log,*)                                           &
           'ERROR in subroutine main_LSLS3: nr_of_LS_coupling.lt.1'
         write(iwrite_log,*)                                           &
          '(not specified a serial number of expansions in LS coupling)'
         write(iwrite_log,*)'program will be terminated'
         stop
      end if
      all_expansions%coupling_expansions(icoupling_nr_in_expansions)%  &
         nr_of_expansions =                                            &
         all_expansions%coupling_expansions(nr_of_LS_coupling)%        &
         nr_of_expansions
      allocate(all_expansions%coupling_expansions(                     &
         icoupling_nr_in_expansions)%expansions(all_expansions%        &
         coupling_expansions(icoupling_nr_in_expansions)%              &
         nr_of_expansions))
      do istate=1, all_expansions%coupling_expansions(                 &
         icoupling_nr_in_expansions)%nr_of_expansions, 1
         all_expansions%coupling_expansions(                           &
           icoupling_nr_in_expansions)%expansions(istate)%nr_of_state= &
           all_expansions%coupling_expansions(nr_of_LS_coupling)%      &
           expansions(istate)%nr_of_state
         call form_cfg_structure(istate)
!        call print_cfg_structure
!     Define the size of mchf expansion in LK coupling:
         nr_of_csfs_LS3=0
         do icfg=1, cfg_structure_LS%nr_of_cfgs
            call form_expansion_cfg_LS(istate, icfg)
!           call print_expansion_cfg_LS
            call count_nr_of_csfs_LS3(cfg_nr_of_csfs_LS3)
            nr_of_csfs_LS3 = nr_of_csfs_LS3 + cfg_nr_of_csfs_LS3
            call delete_expansion_cfg_LS
         end do
         all_expansions%coupling_expansions(                           &
            icoupling_nr_in_expansions)%expansions(istate)%size =      &
            nr_of_csfs_LS3
         all_expansions%coupling_expansions(                           &
            icoupling_nr_in_expansions)%expansions(istate)%            &
            nr_of_state = istate
         allocate(all_expansions%coupling_expansions(                  &
            icoupling_nr_in_expansions)%expansions(istate)%            &
            coeffs(all_expansions%coupling_expansions(                 &
            icoupling_nr_in_expansions)%expansions(istate)%size))
         allocate(all_expansions%coupling_expansions(                  &
            icoupling_nr_in_expansions)%expansions(istate)%            &
            csfs(all_expansions%coupling_expansions(                   &
            icoupling_nr_in_expansions)%expansions(istate)%size))
!
!     Fill up expansions_LK%expansions(istate)
         nr_of_csfs_LS3=0
         sum_of_state = ZERO_dp
         do icfg=1, cfg_structure_LS%nr_of_cfgs
            call form_expansion_cfg_LS(istate, icfg)
            call main_cfg_lsls3(print_level)
            do icsf_LS3=1, expansion_cfg_LS3%size, 1
               nr_of_csfs_LS3 = nr_of_csfs_LS3 + 1
               if(nr_of_csfs_LS3 .gt.                                  &
               all_expansions%coupling_expansions(                     &
               icoupling_nr_in_expansions)%expansions(istate)%size) then
                  write(*,*) 'STOP at subroutine main_lsls3  : ',      &
                             'nr_of_csfs_LS3.gt.',                     &
               'all_expansions%coupling_expansions()%expansions()%size'
                  stop
               end if
!unsing this I am SURE that in "b=a" "b" is "empty" csf ...
               all_expansions%coupling_expansions(                     &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  csfs(nr_of_csfs_LS3)=expansion_cfg_LS3%csfs(icsf_LS3)
               all_expansions%coupling_expansions(                     &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  coeffs(nr_of_csfs_LS3)=expansion_cfg_LS3%            &
                  coeffs(icsf_LS3)
               sum_of_state = sum_of_state + all_expansions%           &
                  coupling_expansions(icoupling_nr_in_expansions)%     &
                  expansions(istate)%coeffs(nr_of_csfs_LS3)*           &
                  all_expansions%coupling_expansions(                  &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  coeffs(nr_of_csfs_LS3)
            end do
            call delete_cfg_expansions
         end do
         if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
            write(*,*)                                                 &
               'POSSIBLE ERROR at subroutine main_LSLS3: ',            &
               'sum_of_state.gt.1: sum_of_state = ', sum_of_state
            write(iwrite_log,*)                                        &
               'POSSIBLE ERROR at subroutine main_LSLK: ',             &
                'sum_of_state.gt.1: sum_of_state = ', sum_of_state
         end if
         if(istate == 1) write(iwrite_log,*)                           &
         'State number   Expansion size   Summation rules'
         write(iwrite_log,*) istate,nr_of_csfs_LS3,'   ',sum_of_state
         call delete_cfg_structure
      end do
!      write(*,*) '  end subroutine main_lsls3'
!
      contains
!
!***********************************************************************
!                                                                      *
         subroutine form_expansion_cfg_LS(istate, icfg)
!                                                                      *
!     This subroutine forms oneconfigurational expansion               *
!     "expansion_cfg_LS" of the state with serial number "istate"      *
!     of configuration with serial number "icfg"                       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
         integer, intent(in) :: istate, icfg
         integer::icsf
         expansion_cfg_LS%size =                                       &
                           cfg_structure_LS%csfs_lists(icfg)%list_size
         expansion_cfg_LS%nr_of_state =                                &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%nr_of_state
         allocate(expansion_cfg_LS%csfs(expansion_cfg_LS%size))
         allocate(expansion_cfg_LS%coeffs(expansion_cfg_LS%size))
         do icsf=1, expansion_cfg_LS%size
            expansion_cfg_LS%coeffs(icsf) =                            &
               all_expansions%coupling_expansions(nr_of_LS_coupling)%  &
               expansions(istate)%coeffs(cfg_structure_LS%             &
               csfs_lists(icfg)%items(icsf))
            expansion_cfg_LS%csfs(icsf) =                              &
               all_expansions%coupling_expansions(nr_of_LS_coupling)%  &
                  expansions(istate)%csfs(cfg_structure_LS%            &
                  csfs_lists(icfg)%items(icsf))
         end do
         end subroutine form_expansion_cfg_LS
!
!***********************************************************************
!                                                                      *
         subroutine delete_expansion_cfg_LS
!                                                                      *
!     This subroutine deallocates the arrays of "expansion_cfg_LS"     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
         integer::error
         integer::icsf
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1, expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg))    &
             deallocate(expansion_cfg_LS%csfs(icsf)%subc_cfg,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%subc))        &
                 deallocate(expansion_cfg_LS%csfs(icsf)%subc,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iM1))         &
                  deallocate(expansion_cfg_LS%csfs(icsf)%iM1,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iM2))         &
                  deallocate(expansion_cfg_LS%csfs(icsf)%iM2,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iJ))         &
                  deallocate(expansion_cfg_LS%csfs(icsf)%iJ,STAT=error)
            end do
            deallocate(expansion_cfg_LS%csfs, STAT=error)
         end if
         deallocate(expansion_cfg_LS%coeffs, STAT=error)
         end subroutine delete_expansion_cfg_LS
!
!subroutine print_cfg_structure
!--------------------------------------------------------------------
! This subroutine prints  "cfg_structure" to the unit "iwrite_log2"
!--------------------------------------------------------------------
!	integer::icfg, icsf
!	write(iwrite_log2,*)'********************************************'
!	write(iwrite_log2,*)'cfg structure: '
!	write(iwrite_log2,*)'istate,  nr_of_cfgs', istate, cfg_structure_LS%nr_of_cfgs
!	!
!	do icfg=1,cfg_structure_LS%nr_of_cfgs
!		write(iwrite_log2,*)'  list_size:', cfg_structure_LS%csfs_lists(icfg)%list_size
!		!write csfs in cfg_structure
!		do icsf=1, cfg_structure_LS%csfs_lists(icfg)%list_size
!			write(iwrite_log2,*)'          csf:', cfg_structure_LS%csfs_lists(icfg)%items(icsf)
!		end do
!	end do
!end subroutine print_cfg_structure
!
!***********************************************************************
!                                                                      *
      end subroutine main_LSLS3
!
      subroutine form_cfg_structure(istate)
!                                                                      *
!     This subroutine forms the structure "cfg_structure_LS"           *
!                                                                      *
!     There assumed that csfs in expansion is ordered so that          *
!     they are in groups of csfs of the same cfg                       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer, intent(in)::istate
      integer::icsf, icfg, inr_of_csfs_of_cfg
      cfg_structure_LS%nr_of_cfgs=cfg_count(istate)
      allocate(cfg_structure_LS%csfs_lists(cfg_structure_LS%           &
         nr_of_cfgs))
      icfg=1
      inr_of_csfs_of_cfg=1
      do icsf=2,all_expansions%coupling_expansions(nr_of_LS_coupling)% &
         expansions(istate)%size
         if(the_same_cfg(all_expansions%coupling_expansions(           &
            nr_of_LS_coupling)%expansions(istate)%csfs(icsf-1),        &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
                                    expansions(istate)%csfs(icsf))) then
            inr_of_csfs_of_cfg = inr_of_csfs_of_cfg + 1
         else
            cfg_structure_LS%csfs_lists(icfg)%list_size =              &
                                                      inr_of_csfs_of_cfg
            allocate(cfg_structure_LS%csfs_lists(icfg)%items(          &
               cfg_structure_LS%csfs_lists(icfg)%list_size))
            inr_of_csfs_of_cfg = 1
            icfg=icfg+1
         end if
      end do
      cfg_structure_LS%csfs_lists(icfg)%list_size=inr_of_csfs_of_cfg
      allocate(cfg_structure_LS%csfs_lists(icfg)%items(                &
         cfg_structure_LS%csfs_lists(icfg)%list_size))
      inr_of_csfs_of_cfg = 1
      do icfg=1, cfg_structure_LS%nr_of_cfgs, 1
         do icsf=1, cfg_structure_LS%csfs_lists(icfg)%list_size, 1
            cfg_structure_LS%csfs_lists(icfg)%items(icsf) =            &
                                                    inr_of_csfs_of_cfg
            inr_of_csfs_of_cfg = inr_of_csfs_of_cfg + 1
         end do
      end do
      end subroutine form_cfg_structure
!
!***********************************************************************
!                                                                      *
      subroutine delete_cfg_structure
!                                                                      *
!     This subroutine deallocates the arrays of "cfg_structure_LS"     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer::error
      integer::icfg
      if(associated(cfg_structure_LS%csfs_lists)) then
         do icfg=1,cfg_structure_LS%nr_of_cfgs
            if(associated(cfg_structure_LS%csfs_lists(icfg)%items))    &
          deallocate(cfg_structure_LS%csfs_lists(icfg)%items,STAT=error)
         end do
         deallocate(cfg_structure_LS%csfs_lists, STAT=error)
      end if
      end subroutine delete_cfg_structure
!
!***********************************************************************
!                                                                      *
      function cfg_count(istate)   result(irez)
!                                                                      *
!     This function counts number of diffferent configurations in      *
!     the expansion of the state "istate" in LS coupling               *
!                                                                      *
!     There assumed that csfs in expansion is ordered so that          *
!     they are in groups of csfs of the same cfg                       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      integer::istate, icsf
      integer::irez
      if(istate.eq.0) then
         irez=0
      else
         if(istate.ne.all_expansions%coupling_expansions(              &
            nr_of_LS_coupling)%expansions(istate)%nr_of_state) then
            write(*,*) 'STOP at subroutine cfg_count module ',         &
               'Coupling_transform_LSLS3: istate.ne.',                 &
               'expansion(nr_of_LS_coupling)(istate)%nr_of_state'
            stop
         end if
         irez=1
         do icsf=2,all_expansions%coupling_expansions(                 &
                             nr_of_LS_coupling)%expansions(istate)%size
            if(.not. the_same_cfg(all_expansions%coupling_expansions(  &
               nr_of_LS_coupling)%expansions(istate)%csfs(icsf-1),     &
               all_expansions%coupling_expansions(nr_of_LS_coupling)%  &
               expansions(istate)%csfs(icsf))) then
               irez=irez+1
            end if
         end do
      end if
      end function
!
      end module Coupling_transform_LSLS3
