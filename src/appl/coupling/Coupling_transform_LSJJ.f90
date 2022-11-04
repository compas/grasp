      module Coupling_transform_LSJJ
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
      use Coupling_transform_cfg_LSJJ
!
      public  :: main_LSJJ
      private :: cfg_count
!
      private :: form_cfg_structure
      private :: delete_cfg_structure
!
      type(cfg_structure), private::cfg_structure_LS
!
!serial number of LS coupling in
!all_axpansions%coupling_expansions()
      integer, private:: nr_of_LS_coupling
!
      contains
!
      subroutine main_LSJJ (icoupling_nr_in_expansions,                &
                                        inr_of_LS_coupling, print_level)
!--------------------------------------------------------------------
! This is managerial subroutine
!--------------------------------------------------------------------
      implicit none
      integer, intent (in):: icoupling_nr_in_expansions
      integer, intent (in):: inr_of_LS_coupling, print_level
      integer::istate
      integer::icfg, icsf_JJ
      real(kind=dp):: sum_of_state
      integer::nr_of_csfs_JJ, cfg_nr_of_csfs_JJ
!
!      write(*,*) '  subroutine main_lsjj'
      write(iwrite_log,*) ''
      write(iwrite_log,*) 'JJ coupling:'
!
      nr_of_LS_coupling = inr_of_LS_coupling
!
      all_expansions%coupling_expansions(icoupling_nr_in_expansions)%  &
         nr_of_expansions = all_expansions%coupling_expansions(        &
         nr_of_LS_coupling)%nr_of_expansions
      allocate(all_expansions%coupling_expansions(                     &
         icoupling_nr_in_expansions)%expansions(                       &
         all_expansions%coupling_expansions(                           &
         icoupling_nr_in_expansions)%nr_of_expansions))
      do istate=1,all_expansions%coupling_expansions(icoupling_nr_in_expansions)%nr_of_expansions,1
!        write(*,*) 'nr_of_cfgs'
!        write(*,'(6x,a7,2x,i2)') 'istate=', istate
        all_expansions%coupling_expansions(icoupling_nr_in_expansions)%&
        expansions(istate)%nr_of_state =                               &
          all_expansions%coupling_expansions(nr_of_LS_coupling)%expansions(istate)%nr_of_state
        call form_cfg_structure(istate)
!call print_cfg_structure
!
!define the size of mchf expansion in JJ coupling:
        nr_of_csfs_JJ=0
        do icfg=1, cfg_structure_LS%nr_of_cfgs
!           write(*,'(9x,a5,2x,i2)') 'icfg=', icfg
           call form_expansion_cfg_LS(istate, icfg)
!call print_expansion_cfg_LS
!
           call count_nr_of_csfs_JJ(cfg_nr_of_csfs_JJ)
           nr_of_csfs_JJ = nr_of_csfs_JJ + cfg_nr_of_csfs_JJ
           call delete_expansion_cfg_LS
        end do !icfg
!
!
        all_expansions%coupling_expansions(icoupling_nr_in_expansions)%expansions(istate)%size=nr_of_csfs_JJ
        all_expansions%coupling_expansions(icoupling_nr_in_expansions)%expansions(istate)%nr_of_state=istate
!
!allocate all_expansions%coupling_expansions(icoupling_nr_in_expansions)%expansions(istate)
        allocate(                                                      &
       all_expansions%coupling_expansions(icoupling_nr_in_expansions)% &
        expansions(istate)%coeffs(all_expansions%coupling_expansions(  &
        icoupling_nr_in_expansions)%expansions(istate)%size))
        allocate(                                                      &
       all_expansions%coupling_expansions(icoupling_nr_in_expansions)% &
        expansions(istate)%csfs(all_expansions%coupling_expansions(    &
        icoupling_nr_in_expansions)%expansions(istate)%size))
!
!fill up expansions_JJ%expansions(istate)
        nr_of_csfs_JJ=0
        sum_of_state = 0.
        do icfg=1, cfg_structure_LS%nr_of_cfgs
!write(iwrite_log,'(9x,a5,2x,i2)') 'icfg=', icfg
           call form_expansion_cfg_LS(istate, icfg)
           call main_cfg_lsjj(print_level)
           do icsf_JJ=1, expansion_cfg_JJ%size, 1
              nr_of_csfs_JJ = nr_of_csfs_JJ + 1
              if(nr_of_csfs_JJ .gt.                                    &
                      all_expansions%coupling_expansions(              &
               icoupling_nr_in_expansions)%expansions(istate)%size) then
                 write(*,*) 'STOP at subroutine main_lsjj: ',          &
                 'nr_of_csfs_JJ.gt.all_expansions%coupling_expansion', &
                 '()%expansions()%size'
                 stop
              end if
!unsing this I am SURE that in "b=a" "b" is "empty" csf ...
              all_expansions%coupling_expansions(                      &
                 icoupling_nr_in_expansions)%expansions(istate)%csfs(  &
                 nr_of_csfs_JJ)=expansion_cfg_JJ%csfs(icsf_JJ)
              all_expansions%coupling_expansions(                      &
                 icoupling_nr_in_expansions)%expansions(istate)%       &
                 coeffs(nr_of_csfs_JJ)=expansion_cfg_JJ%coeffs(icsf_JJ)
              sum_of_state = sum_of_state +                            &
                 all_expansions%coupling_expansions(                   &
                 icoupling_nr_in_expansions)%expansions(istate)%       &
                 coeffs(nr_of_csfs_JJ)* &
                 all_expansions%coupling_expansions(                   &
                 icoupling_nr_in_expansions)%expansions(istate)%       &
                 coeffs(nr_of_csfs_JJ)
           end do
           call delete_cfg_expansions
        end do
        if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
           write(*,*) 'POSSIBLE ERROR at subroutine main_LSJJ: ',      &
              'sum_of_state.gt.1: sum_of_state = ', sum_of_state
           write(*,*) 'State number = ', istate
           write(iwrite_log,*) 'POSSIBLE ERROR at subroutine ',        &
              'main_LSJJ: sum_of_state.gt.1: sum_of_state = ',         &
              sum_of_state
           write(iwrite_log,*) 'State number = ', istate
           stop
        end if
        if( istate == 1) write(iwrite_log,*)                           &
        'State number   Expansion size   Summation rules'
        write(iwrite_log,*) istate,nr_of_csfs_JJ,'   ',sum_of_state
        call delete_cfg_structure
      end do
!      write(*,*) '  end subroutine main_lsjj'
!
      contains
!
         subroutine form_expansion_cfg_LS(istate, icfg)
!--------------------------------------------------------------------
! This subroutine forms oneconfigurational
! expansion "expansion_cfg_LS"
! of the state with serial number "istate"
! of configuration with serial number "icfg"
!--------------------------------------------------------------------
         integer, intent(in) :: istate, icfg
         integer::icsf
         expansion_cfg_LS%size=cfg_structure_LS%csfs_lists(icfg)%list_size
         expansion_cfg_LS%nr_of_state = all_expansions%coupling_expansions(nr_of_LS_coupling)%expansions(istate)%nr_of_state
         allocate(expansion_cfg_LS%csfs(expansion_cfg_LS%size))
         allocate(expansion_cfg_LS%coeffs(expansion_cfg_LS%size))
         do icsf=1, expansion_cfg_LS%size
            expansion_cfg_LS%coeffs(icsf) =                            &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%coeffs(cfg_structure_LS%                &
            csfs_lists(icfg)%items(icsf))
            expansion_cfg_LS%csfs(icsf) =                              &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%csfs(cfg_structure_LS%                  &
            csfs_lists(icfg)%items(icsf))
         end do
         end subroutine form_expansion_cfg_LS
!
         subroutine delete_expansion_cfg_LS
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of "expansion_cfg_LS"
!--------------------------------------------------------------------
         integer::error
         integer::icsf
!
!      write(*,*) '  delete_expansion_cfg_LS'
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1, expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg))   &
            deallocate(expansion_cfg_LS%csfs(icsf)%subc_cfg,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%subc))       &
            deallocate(expansion_cfg_LS%csfs(icsf)%subc,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iM1))        &
            deallocate(expansion_cfg_LS%csfs(icsf)%iM1,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iM2))        &
            deallocate(expansion_cfg_LS%csfs(icsf)%iM2,STAT=error)
               if(associated(expansion_cfg_LS%csfs(icsf)%iJ))         &
            deallocate(expansion_cfg_LS%csfs(icsf)%iJ, STAT=error)
            end do
            deallocate(expansion_cfg_LS%csfs, STAT=error)
         end if !assoctated(expansion_cfg_LS%csfs)
!
         deallocate(expansion_cfg_LS%coeffs, STAT=error)
!      write(*,*) '  delete_expansion_cfg_LS'
         end subroutine delete_expansion_cfg_LS
!
!subroutine print_expansion_cfg_LS
!--------------------------------------------------------------------
! This subroutine
!--------------------------------------------------------------------
!	integer::icsf, isubc
!	integer::il, in, iN_big, inr, Li, Si, iM1, iM2
!	write(iwrite_log2,*)'  -------- cfg LS expansion ----------'
!	write(iwrite_log2,*)'  size:',expansion_cfg_LS%size
!	if(associated(expansion_cfg_LS%csfs)) then
!	do icsf=1, expansion_cfg_LS%size, 1
!		if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
!		do isubc=1,expansion_cfg_LS%csfs(icsf)%nosubc,1
!			iN_big= expansion_cfg_LS%csfs(icsf)%subc_cfg(isubc)%iN_big
!			il = expansion_cfg_LS%csfs(icsf)%subc_cfg(isubc)%il
!			in = expansion_cfg_LS%csfs(icsf)%subc_cfg(isubc)%in
!			write(iwrite_log2,'(9x,i2,3x,i2,3x,i2,3x,i2)') isubc,in,il,iN_big
!		end do
!		else
!		   write(iwrite_log2,'(9x,a51)')'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
!		end if !	associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)
!		write(iwrite_log2,'(8x,a5)')'state'
!		write(iwrite_log2,'(9x,a29)')'subc  nr   L    S   L_i   S_i'
!		if(associated(expansion_cfg_LS%csfs(icsf)%subc).and.associated(expansion_cfg_LS%csfs(icsf)%iM1).and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
!		do isubc=1,expansion_cfg_LS%csfs(icsf)%nosubc,1
!			inr= expansion_cfg_LS%csfs(icsf)%subc(isubc)%inr
!			Li= expansion_cfg_LS%csfs(icsf)%subc(isubc)%iL
!			Si= expansion_cfg_LS%csfs(icsf)%subc(isubc)%iS
!			iM1= expansion_cfg_LS%csfs(icsf)%iM1(isubc)
!			iM2= expansion_cfg_LS%csfs(icsf)%iM2(isubc)
!			write(iwrite_log2,'(9x,i2,3x,i2,3x,i2,3x,i2,3x,i2,3x,i2)')isubc,inr,Li,Si,iM1,iM2
!		end do
!		else
!		   write(iwrite_log2,'(9x,a63)')'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
!		end if !	associated(expansion_cfg_LS%csfs(icsf)%subc).and.associated(expansion_cfg_JJ%csfs(icsf)%iM1).and.associated(expansion_cfg_JJ%csfs(icsf)%iM2)) then
!
!		if(associated(expansion_cfg_LS%coeffs)) then
!			write(iwrite_log2,'(4x,f9.7)') expansion_cfg_LS%coeffs(icsf)
!		else
!			write(iwrite_log2,'(4x,a33)')'expansion_cfg_LS%coeffs NOT associated'
!		end if !associated(expansion_cfg_LS%coeffs)
!	end do
!	else
!		write(iwrite_log2,*)'expansion_cfg_LS%csfs NOT associated'
!	end if
!	!
!	write(iwrite_log2,*)'  ------ END cfg LS expansion --------'
!end subroutine print_expansion_cfg_LS
!
!subroutine print_cfg_structure
!	integer::icfg, icsf
!	write(iwrite_log2,*)'********************************************'
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
      end subroutine main_LSJJ
!
!---------  subroutine form_cfg_structure  ------
!
      subroutine form_cfg_structure(istate)
!--------------------------------------------------------------------
! This subroutine forms the structure "cfg_structure_LS"
!
!there assumed that csfs in expansion is ordered so that
!they are in groups of csfs of the same cfg
!--------------------------------------------------------------------
      implicit none
      integer, intent(in)::istate
      integer::icsf, icfg, inr_of_csfs_of_cfg
      cfg_structure_LS%nr_of_cfgs=cfg_count(istate)
      allocate(cfg_structure_LS%csfs_lists(cfg_structure_LS%nr_of_cfgs))
      icfg=1
      inr_of_csfs_of_cfg=1
      do icsf=2,all_expansions%coupling_expansions(nr_of_LS_coupling)%expansions(istate)%size
         if(the_same_cfg(all_expansions%coupling_expansions(           &
            nr_of_LS_coupling)%expansions(istate)%csfs(icsf-1),        &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%csfs(icsf))) then
            inr_of_csfs_of_cfg = inr_of_csfs_of_cfg + 1
         else
            cfg_structure_LS%csfs_lists(icfg)%list_size=inr_of_csfs_of_cfg
            allocate(cfg_structure_LS%csfs_lists(icfg)%items(cfg_structure_LS%csfs_lists(icfg)%list_size))
            inr_of_csfs_of_cfg = 1
            icfg=icfg+1
         end if
      end do
      cfg_structure_LS%csfs_lists(icfg)%list_size=inr_of_csfs_of_cfg
      allocate(cfg_structure_LS%csfs_lists(icfg)%items(cfg_structure_LS%csfs_lists(icfg)%list_size))
      inr_of_csfs_of_cfg = 1
      do icfg=1, cfg_structure_LS%nr_of_cfgs, 1
         do icsf=1, cfg_structure_LS%csfs_lists(icfg)%list_size, 1
            cfg_structure_LS%csfs_lists(icfg)%items(icsf)=inr_of_csfs_of_cfg
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
! This subroutine deallocates the arrays of "cfg_structure_LS"
!--------------------------------------------------------------------
      implicit none
      integer::error
      integer::icfg
      if(associated(cfg_structure_LS%csfs_lists)) then
         do icfg=1,cfg_structure_LS%nr_of_cfgs
            if(associated(cfg_structure_LS%csfs_lists(icfg)%items)) deallocate(cfg_structure_LS%csfs_lists(icfg)%items, STAT=error)
         end do
         deallocate(cfg_structure_LS%csfs_lists, STAT=error)
      end if
!
      end subroutine delete_cfg_structure
!
!--------  function cfg_count    ----------------
!
      function cfg_count(istate)   result(irez)
!--------------------------------------------------------------------
! This function counts number of diffferent
! configurations in the expansion of the state "istate"
! in LS coupling
!
!there assumed that csfs in expansion is ordered so that
!they are in groups of csfs of the same cfg
!--------------------------------------------------------------------
      implicit none
      integer::istate, icsf
      integer::irez
      if(istate.eq.0) then
         irez=0
      else
         if(istate .ne.                                                &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%nr_of_state) then
            write(*,*) 'STOP at subroutine cfg_count module ',         &
            'Coupling_transform_LSJJ: ',                               &
            'istate.ne.expansion()()%nr_of_state'
            write(*,*) istate,                                         &
            all_expansions%coupling_expansions(nr_of_LS_coupling)%     &
            expansions(istate)%nr_of_state
            stop
         end if
         irez=1
         do icsf=2,all_expansions%coupling_expansions(                 &
                             nr_of_LS_coupling)%expansions(istate)%size
            if(.not.                                                   &
               the_same_cfg(all_expansions%coupling_expansions(        &
               nr_of_LS_coupling)%expansions(istate)%csfs(icsf-1),     &
               all_expansions%coupling_expansions(nr_of_LS_coupling)%  &
                                   expansions(istate)%csfs(icsf))) then
               irez=irez+1
            end if
         end do
      end if
      end function
!
end module Coupling_transform_LSJJ
