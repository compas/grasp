      module Coupling_transform_JJJK
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
      use Coupling_transform_cfg_JJJK
!
      public  :: main_JJJK
      private :: cfg_count
!
      private :: form_cfg_structure
      private :: delete_cfg_structure
!
      type(cfg_structure), private::cfg_structure_JJ
!
      integer, private :: nr_of_JJ_coupling
!
      contains
!
      subroutine main_JJJK (icoupling_nr_in_expansions,                &
                                        iinr_of_JJ_coupling,print_level)
!--------------------------------------------------------------------
! This is managerial subroutine
!--------------------------------------------------------------------
      implicit none
      integer, intent (in):: icoupling_nr_in_expansions
      integer, intent (in):: iinr_of_JJ_coupling, print_level
      integer::istate
      integer::icfg, icsf_JK
      integer::nr_of_csfs_JK, cfg_nr_of_csfs_JK
      real(kind=dp)::sum_of_state
!      write(*,*) '  subroutine main_JJJK'
      write(iwrite_log,*) ''
      write(iwrite_log,*) 'JK coupling:'
      nr_of_JJ_coupling = iinr_of_JJ_coupling
      if(nr_of_JJ_coupling.lt.1) then
         write(*,*)'ERROR in subroutine main_JJJK: ',                  &
            'nr_of_JJ_coupling.lt.1'
         write(*,*)'(not specified a serial number of expansions in ', &
            'JJ coupling); nr_of_JJ_coupling=', nr_of_JJ_coupling
         write(*,*)'program will be terminated'
         write(iwrite_log,*)'ERROR in subroutine main_JJJK: ',         &
            'nr_of_JJ_coupling.lt.1'
         write(iwrite_log,*)'(not specified a serial number of ',      &
            'expansions in JJ coupling); nr_of_JJ_coupling=',          &
            nr_of_JJ_coupling
         write(iwrite_log,*)'program will be terminated'
         stop
      end if
      if(coupling_descriptions(all_expansions%coupling_expansions(     &
                 nr_of_JJ_coupling)%icoupling)%short_name.ne.'JJ') then
         write(*,*) 'ERROR AT subroutine main_JJJK: ',                 &
           'coupling_descriptions(all_expansions%coupling_expansions(',&
           'inr_of_JJ_coupling)%icoupling)%short_name.ne.JJ'
         write(*,*) '(possibly - using JJ-JK transformation without ', &
           'JJ coupling, or coupling_descriptions(...) not correctly ',&
           'filled...)'
         write(*,*) 'program will be terminated'
         write(iwrite_log,*) 'ERROR AT subroutine main_JJJK: ',        &
            'coupling_descriptions(all_expansions%',                   &
            'coupling_expansions()%icoupling)%short_name.ne.JJ'
         write(iwrite_log,*) '(possibly - using JJ-JK transformation ',&
            'without JJ coupling, or coupling_descriptions(...) ',     &
            'not correctly filled...)'
         write(iwrite_log,*) 'program will be terminated'
         stop
      end if
      all_expansions%coupling_expansions(icoupling_nr_in_expansions)%  &
         nr_of_expansions = all_expansions%coupling_expansions(        &
         nr_of_JJ_coupling)%nr_of_expansions
      allocate(all_expansions%coupling_expansions(                     &
         icoupling_nr_in_expansions)%expansions(all_expansions%        &
         coupling_expansions(icoupling_nr_in_expansions)%              &
         nr_of_expansions))
      do istate=1, all_expansions%coupling_expansions(                 &
         icoupling_nr_in_expansions)%nr_of_expansions, 1
         all_expansions%coupling_expansions(                           &
            icoupling_nr_in_expansions)%expansions(istate)%nr_of_state &
            = all_expansions%coupling_expansions(nr_of_JJ_coupling)%   &
            expansions(istate)%nr_of_state
         call form_cfg_structure(istate)
         nr_of_csfs_JK=0
         do icfg=1, cfg_structure_JJ%nr_of_cfgs
            call form_expansion_cfg_JJ(istate, icfg)
            call count_nr_of_csfs_JK(cfg_nr_of_csfs_JK)
            nr_of_csfs_JK = nr_of_csfs_JK + cfg_nr_of_csfs_JK
            call delete_expansion_cfg_JJ
         end do
         all_expansions%coupling_expansions(                           &
            icoupling_nr_in_expansions)%expansions(istate)%size =      &
            nr_of_csfs_JK
         all_expansions%coupling_expansions(                           &
            icoupling_nr_in_expansions)%expansions(istate)%            &
            nr_of_state = istate
         allocate(all_expansions%coupling_expansions(                  &
            icoupling_nr_in_expansions)%expansions(istate)%coeffs(     &
            all_expansions%coupling_expansions(                        &
            icoupling_nr_in_expansions)%expansions(istate)%size))
         allocate(                                                     &
        all_expansions%coupling_expansions(icoupling_nr_in_expansions)%&
         expansions(istate)%csfs(all_expansions%coupling_expansions(   &
         icoupling_nr_in_expansions)%expansions(istate)%size))
         nr_of_csfs_JK=0
         sum_of_state = ZERO_dp
         do icfg=1, cfg_structure_JJ%nr_of_cfgs
            call form_expansion_cfg_JJ(istate, icfg)
            call main_cfg_JJJK(print_level)
            do icsf_JK=1, expansion_cfg_JK%size, 1
               nr_of_csfs_JK = nr_of_csfs_JK + 1
               if(nr_of_csfs_JK.gt.all_expansions%coupling_expansions( &
                   icoupling_nr_in_expansions)%expansions(istate)%size)&
                                                             then
                   write(*,*) 'STOP at subroutine main_JJjj  : ',      &
                    'nr_of_csfs_JK.gt.',                               &
               'all_expansions%coupling_expansions()%expansions()%size'
                   stop
               end if
               all_expansions%coupling_expansions(                     &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  csfs(nr_of_csfs_JK)=expansion_cfg_JK%csfs(icsf_JK)
               all_expansions%coupling_expansions(                     &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  coeffs(nr_of_csfs_JK)=expansion_cfg_JK%coeffs(icsf_JK)
               sum_of_state = sum_of_state +                           &
                  all_expansions%coupling_expansions(                  &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  coeffs(nr_of_csfs_JK) *                              &
                  all_expansions%coupling_expansions(                  &
                  icoupling_nr_in_expansions)%expansions(istate)%      &
                  coeffs(nr_of_csfs_JK)
            end do
            call delete_cfg_expansions
         end do
         if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
            write(*,*) 'POSSIBLE ERROR at subroutine main_JJJK: ',     &
               'sum_of_state.gt.1: sum_of_state = ', sum_of_state
            write(iwrite_log,*)                                        &
               'POSSIBLE ERROR at subroutine main_JJJK: ',             &
               'sum_of_state.gt.1: sum_of_state = ', sum_of_state
         end if
         if(istate == 1) write(iwrite_log,*)                           &
         'State number   Expansion size   Summation rules'
         write(iwrite_log,*) istate,nr_of_csfs_JK,'   ',sum_of_state
         call delete_cfg_structure
      end do
!      write(*,*) '  end subroutine main_JJJK'
!
      contains
!
         subroutine form_expansion_cfg_JJ(istate, icfg)
!--------------------------------------------------------------------
! This subroutine forms oneconfigurational
! expansion "expansion_cfg_JJ"
! of the state with serial number "istate"
! of configuration with serial number "icfg"
!--------------------------------------------------------------------
         integer, intent(in) :: istate, icfg
         integer::icsf
         expansion_cfg_JJ%size =                                       &
                           cfg_structure_JJ%csfs_lists(icfg)%list_size
         expansion_cfg_JJ%nr_of_state =                                &
            all_expansions%coupling_expansions(nr_of_JJ_coupling)%     &
            expansions(istate)%nr_of_state
         allocate(expansion_cfg_JJ%csfs(expansion_cfg_JJ%size))
         allocate(expansion_cfg_JJ%coeffs(expansion_cfg_JJ%size))
         do icsf=1, expansion_cfg_JJ%size
            expansion_cfg_JJ%coeffs(icsf) =                            &
               all_expansions%coupling_expansions(nr_of_JJ_coupling)%  &
               expansions(istate)%coeffs(cfg_structure_JJ%csfs_lists(  &
               icfg)%items(icsf))
            expansion_cfg_JJ%csfs(icsf) =                              &
               all_expansions%coupling_expansions(nr_of_JJ_coupling)%  &
               expansions(istate)%csfs(cfg_structure_JJ%csfs_lists(    &
               icfg)%items(icsf))
         end do
         end subroutine form_expansion_cfg_JJ
!
         subroutine delete_expansion_cfg_JJ
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of "expansion_cfg_JJ"
!--------------------------------------------------------------------
         integer::error
         integer::icsf
         if(associated(expansion_cfg_JJ%csfs)) then
            do icsf=1, expansion_cfg_JJ%size
               if(associated(expansion_cfg_JJ%csfs(icsf)%subc_cfg))    &
             deallocate(expansion_cfg_JJ%csfs(icsf)%subc_cfg,STAT=error)
               if(associated(expansion_cfg_JJ%csfs(icsf)%subc))        &
                 deallocate(expansion_cfg_JJ%csfs(icsf)%subc,STAT=error)
               if(associated(expansion_cfg_JJ%csfs(icsf)%iM1))         &
                 deallocate(expansion_cfg_JJ%csfs(icsf)%iM1, STAT=error)
               if(associated(expansion_cfg_JJ%csfs(icsf)%iM2))         &
                 deallocate(expansion_cfg_JJ%csfs(icsf)%iM2, STAT=error)
               if(associated(expansion_cfg_JJ%csfs(icsf)%iJ))          &
                 deallocate(expansion_cfg_JJ%csfs(icsf)%iJ, STAT=error)
            end do
            deallocate(expansion_cfg_JJ%csfs, STAT=error)
         end if
         deallocate(expansion_cfg_JJ%coeffs, STAT=error)
         end subroutine delete_expansion_cfg_JJ
!
!subroutine print_cfg_structure
!--------------------------------------------------------------------
! This subroutine prints the "cfg_structure_JJ" to the unit "iwrite_log2"
!--------------------------------------------------------------------
!	integer::icfg, icsf
!	write(iwrite_log2,*)'********************************************'
!	write(iwrite_log2,*)'cfg structure: '
!	write(iwrite_log2,*)'istate,  nr_of_cfgs', istate, cfg_structure_JJ%nr_of_cfgs
!	!
!	do icfg=1,cfg_structure_JJ%nr_of_cfgs
!		write(iwrite_log2,*)'  list_size:', cfg_structure_JJ%csfs_lists(icfg)%list_size
!		!write csfs in cfg_structure
!		do icsf=1, cfg_structure_JJ%csfs_lists(icfg)%list_size
!			write(iwrite_log2,*)'          csf:', cfg_structure_JJ%csfs_lists(icfg)%items(icsf)
!		end do
!	end do
!end subroutine print_cfg_structure
!
      end subroutine main_JJJK

!---------  subroutine form_cfg_structure  ------
!
      subroutine form_cfg_structure(istate)
!--------------------------------------------------------------------
! This subroutine forms the structure "cfg_structure_JJ"
!--------------------------------------------------------------------
      implicit none
      integer, intent(in)::istate
      integer::icsf, icfg, inr_of_csfs_of_cfg
      cfg_structure_JJ%nr_of_cfgs=cfg_count(istate)
      allocate(cfg_structure_JJ%csfs_lists(cfg_structure_JJ%           &
                                                           nr_of_cfgs))
      icfg=1
      inr_of_csfs_of_cfg=1
      do icsf=2,all_expansions%coupling_expansions(nr_of_JJ_coupling)% &
                                                expansions(istate)%size
         if(the_same_cfg(all_expansions%coupling_expansions(           &
            nr_of_JJ_coupling)%expansions(istate)%csfs(icsf-1),        &
            all_expansions%coupling_expansions(nr_of_JJ_coupling)%     &
            expansions(istate)%csfs(icsf))) then
            inr_of_csfs_of_cfg = inr_of_csfs_of_cfg + 1
         else
            cfg_structure_JJ%csfs_lists(icfg)%list_size =              &
               inr_of_csfs_of_cfg
            allocate(cfg_structure_JJ%csfs_lists(icfg)%items(          &
              cfg_structure_JJ%csfs_lists(icfg)%list_size))
            inr_of_csfs_of_cfg = 1
            icfg=icfg+1
         end if
      end do
      cfg_structure_JJ%csfs_lists(icfg)%list_size=inr_of_csfs_of_cfg
      allocate(cfg_structure_JJ%csfs_lists(icfg)%items(                &
        cfg_structure_JJ%csfs_lists(icfg)%list_size))
      inr_of_csfs_of_cfg = 1
      do icfg=1, cfg_structure_JJ%nr_of_cfgs, 1
         do icsf=1, cfg_structure_JJ%csfs_lists(icfg)%list_size, 1
            cfg_structure_JJ%csfs_lists(icfg)%items(icsf) =            &
               inr_of_csfs_of_cfg
            inr_of_csfs_of_cfg = inr_of_csfs_of_cfg + 1
         end do
      end do
!
      end subroutine form_cfg_structure
!
!---------  subroutine delete_cfg_structure  ------
!
      subroutine delete_cfg_structure
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of the structure "cfg_structure_JJ"
!--------------------------------------------------------------------
      implicit none
      integer::error
      integer::icfg
      if(associated(cfg_structure_JJ%csfs_lists)) then
         do icfg=1,cfg_structure_JJ%nr_of_cfgs
            if(associated(cfg_structure_JJ%csfs_lists(icfg)%items))    &
               deallocate(cfg_structure_JJ%csfs_lists(icfg)%items,     &
                                                            STAT=error)
         end do
         deallocate(cfg_structure_JJ%csfs_lists, STAT=error)
      end if
      end subroutine delete_cfg_structure
!
!--------  function cfg_count    ----------------
!
      function cfg_count(istate)   result(irez)
!--------------------------------------------------------------------
! This function counts number of diffferent
! configurations in the expansion of the state "istate"
!--------------------------------------------------------------------
      implicit none
      integer::istate, icsf
      integer::irez
      if(istate.eq.0) then
         irez=0
      else
         if(istate.ne.all_expansions%coupling_expansions(              &
            nr_of_JJ_coupling)%expansions(istate)%nr_of_state) then
            write(*,*)                                                 &
        'STOP at subroutine cfg_count module Coupling_transform_JJJK:',&
           ' istate.ne.expansion(nr_of_JJ_coupling)(istate)%nr_of_state'
            stop
         end if
         irez=1
         do icsf=2,all_expansions%coupling_expansions(                 &
                             nr_of_JJ_coupling)%expansions(istate)%size
!there assumed that csfs in expansion is ordered so that
!they are in groups of csfs of teh same cfg
            if(.not. the_same_cfg(all_expansions%coupling_expansions(  &
               nr_of_JJ_coupling)%expansions(istate)%csfs(icsf-1),     &
               all_expansions%coupling_expansions(nr_of_JJ_coupling)%  &
               expansions(istate)%csfs(icsf))) then
               irez=irez+1
            end if
         end do
      end if
      end function
!
      end module
