      module Coupling_transform_cfg_LSJJ
!
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
!
      public :: main_cfg_lsjj
      public :: count_nr_of_csfs_JJ
      public :: delete_cfg_expansions
      private:: form_list_nomach_csfs
      private:: form_csfs_JJ
      private:: matrix_LS_JJ
      private:: dotransform_lsjj
      private:: print_cfg_lsjj
!
      type::Ji_lists
!integer::nr_of_subc
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
      type(expansion), public :: expansion_cfg_JJ
!
      type(list), private ::nomach_csfs_ls
      type(cfg_Ji_lists), private::cfg_Ji_structure
!
      contains
!
      subroutine main_cfg_lsjj(print_level)
!subroutine main_cfg_lsjj(in_expansion_LS, expansion_JJ)
!--------------------------------------------------------------------
! This is managerial subroutine
!--------------------------------------------------------------------
      implicit none
      integer, intent(in) :: print_level
      integer::error
      integer::itype
!
!write(*,*) '    subroutine main_cfg_lsjj'
!!write(iwrite_log,*) '    subroutine main_cfg_lsjj'
!
      expansion_cfg_JJ%nr_of_state=expansion_cfg_LS%nr_of_state
!
      call form_list_nomach_csfs
!
      itype=0
      call form_csfs_JJ(itype)
!
      call dotransform_lsjj
!
      if(print_level.gt.1) call print_cfg_lsjj(2) ! 2-means print LS and JJ expansions
!
      deallocate(nomach_csfs_ls%items, STAT = error)
!
!call delete_cfg_expansions
!
!write(*,*) '    subroutine main_cfg_lsjj'
!write(iwrite_log,*) '    end subroutine main_cfg_lsjj'
      end subroutine main_cfg_lsjj
!
      subroutine count_nr_of_csfs_JJ(nr_of_csfs)
!--------------------------------------------------------------------
! This subroutine counts the number of oneconfigurational
! expansion in JJ coupling
! (corresponding to the one in the JJ coupling "expansion_cfg_LS")
!--------------------------------------------------------------------
      implicit none
      integer, intent(out):: nr_of_csfs
      integer::error
      integer::itype
!
!write(*,*) '    subroutine main_cfg_lsjj'
!write(iwrite_log,*) '    count_nr_of_csfs_JJ'
!call print_cfg_lsjj(1)
!
      expansion_cfg_JJ%nr_of_state=expansion_cfg_LS%nr_of_state
!
      call form_list_nomach_csfs
      itype=1
      call form_csfs_JJ(itype)
!
      nr_of_csfs = expansion_cfg_JJ%size
!
      if(associated(nomach_csfs_ls%items))              &
               deallocate(nomach_csfs_ls%items,STAT=error)
!
!write(*,*) '    subroutine main_cfg_lsjj'
!write(iwrite_log,*) '    end count_nr_of_csfs_JJ'
      end subroutine count_nr_of_csfs_JJ
!
      subroutine delete_cfg_expansions
!--------------------------------------------------------------------
! This subroutine deallocates the arrays of
! oneconfigurational expansions "expansion_cfg_JK" and
! "expansion_cfg_LS"
!--------------------------------------------------------------------
      integer::icsf
!
!write(iwrite_log,'(6x,a21)') 'delete_cfg_expansions'
!
!delete expansion_cfg_JJ
      if(associated(expansion_cfg_JJ%coeffs))                        &
                         deallocate(expansion_cfg_JJ%csfs)
      if(associated(expansion_cfg_JJ%csfs)) then
         do icsf=1, expansion_cfg_JJ%size
            if(associated(expansion_cfg_JJ%csfs(icsf)%subc_cfg))     &
                     deallocate(expansion_cfg_JJ%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_JJ%csfs(icsf)%subc))         &
                     deallocate(expansion_cfg_JJ%csfs(icsf)%subc)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iM1))          &
                     deallocate(expansion_cfg_JJ%csfs(icsf)%iM1)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iM2))          &
                     deallocate(expansion_cfg_JJ%csfs(icsf)%iM2)
            if(associated(expansion_cfg_JJ%csfs(icsf)%iJ))          &
                     deallocate(expansion_cfg_JJ%csfs(icsf)%iJ)
         end do
         deallocate(expansion_cfg_JJ%csfs)
      end if
      expansion_cfg_JJ%size=0
      expansion_cfg_JJ%nr_of_state=0
!
!delete expansion_cfg_LS
      if(associated(expansion_cfg_LS%coeffs))                        &
                     deallocate(expansion_cfg_LS%csfs)
      if(associated(expansion_cfg_LS%csfs)) then
         do icsf=1, expansion_cfg_LS%size
            if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg))     &
                     deallocate(expansion_cfg_LS%csfs(icsf)%subc_cfg)
            if(associated(expansion_cfg_LS%csfs(icsf)%subc))         &
                     deallocate(expansion_cfg_LS%csfs(icsf)%subc)
            if(associated(expansion_cfg_LS%csfs(icsf)%iM1))          &
                     deallocate(expansion_cfg_LS%csfs(icsf)%iM1)
            if(associated(expansion_cfg_LS%csfs(icsf)%iM2))          &
                     deallocate(expansion_cfg_LS%csfs(icsf)%iM2)
            if(associated(expansion_cfg_LS%csfs(icsf)%iJ))          &
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
!-------  subroutine form_list_nomach_csfs   ------------
!
      subroutine form_list_nomach_csfs
!--------------------------------------------------------------------
! This subroutine form the list of serial numbers
! (in "expansion_cfg_LS") of "nonequivalent"
! LS coupling csfs
!(for notion of the "equivalency" see the
! description in the program)
!--------------------------------------------------------------------
      implicit none
      integer:: inomach_counter !,ipast_csf_nr
      integer::icsf_LS, inew
      integer, dimension(:), pointer :: temp_list
      logical:: new_csf
!
      allocate(temp_list(expansion_cfg_LS%size))
!
      inomach_counter=0
!ipast_csf_nr = 0
      do icsf_LS=1,expansion_cfg_LS%size,1
         new_csf = .true.
         do inew = 1,inomach_counter,1
            if(the_same_subc(expansion_cfg_LS%csfs(icsf_LS),           &
              expansion_cfg_LS%csfs(temp_list(inew)))) new_csf = .false.
         end do
         if(new_csf) then
            inomach_counter = inomach_counter + 1
            if(inomach_counter.gt.expansion_cfg_LS%size) then
               write(*,*)'ERROR at subroutine form_list_nomach_csfs: ',&
                         'inomach_counter.gt.expansion_cfg_LS%size'
               stop
            end if
            temp_list(inomach_counter) = icsf_LS
         end if
      end do !icsf_LS
!
      nomach_csfs_LS%list_size = inomach_counter
      allocate(nomach_csfs_LS%items(nomach_csfs_LS%list_size))
      inomach_counter=0
      do inomach_counter=1, nomach_csfs_LS%list_size, 1
         nomach_csfs_LS%items(inomach_counter) =       &
                              temp_list(inomach_counter)
      end do !icsf_LS
      deallocate(temp_list)
      end subroutine form_list_nomach_csfs
!
!***********************************************************************
!                                                                      *
      subroutine form_csfs_JJ(itype)
!                                                                      *
!     This subroutine forms the oneconfigurational                     *
!     expansion in JJ coupling "expansion_cfg_JJ",                     *
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
                      !of csfs_JJ only,
                      !itype.ne.1 - perform full
                      !calculation
      integer,dimension(:),pointer::Ji
      integer,dimension(:),pointer::J_i
      integer:: inr_of_csfs_JJ
      integer:: icsf_LS, icsf_JJ
      integer:: i1, i2
      integer:: J12_min, J12_max, J12
      integer:: isubc, icsf
!integer:: iLi, iSi, inri, inui
      integer::iJ_total
      integer::error
      integer::inr_of_csf_in_expansion_LS
!
!write(*,*) '      subroutine form_csfs_JJ'
!write(iwrite_log,*) '      subroutine form_csfs_JJ'
!
      iJ_total=states%states(expansion_cfg_LS%nr_of_state)%J
!
      call form_Ji_structure
!call show_Ji_structure(55)
      call define_number_of_csfs_JJ(inr_of_csfs_JJ)
!
      expansion_cfg_JJ%size=inr_of_csfs_JJ
      if(itype.ne.1) then
         allocate(expansion_cfg_JJ%csfs(expansion_cfg_JJ%size))
         allocate(expansion_cfg_JJ%coeffs(expansion_cfg_JJ%size))
!
         allocate(J_i(cfg_Ji_structure%nr_of_subc))
         allocate(Ji(cfg_Ji_structure%nr_of_subc))
!
         icsf_JJ=0
         do icsf_LS=1,cfg_Ji_structure%nr_of_nomach_csfs
!define csf number in expansion_cfg_LS
            inr_of_csf_in_expansion_LS=nomach_csfs_ls%items(icsf_LS)
            do i1=1,cfg_Ji_structure%csfs_Ji(icsf_LS)%Ji(1)%list_size
!write(iwrite_lsjj,*)'       i1  ',i1
               Ji(1)=cfg_Ji_structure%csfs_Ji(icsf_LS)%Ji(1)%items(i1)
               J_i(1)=Ji(1)  !tai teisinga tik subc=1 atveju
!write(iwrite_lsjj,*)'       J1  ',J1
               if(cfg_Ji_structure%nr_of_subc.gt.1) then
                 do i2=1,cfg_Ji_structure%csfs_Ji(icsf_LS)%Ji(2)%      &
                                                              list_size
!write(iwrite_lsjj,*)'          i2  ',i2
                   Ji(2)=cfg_Ji_structure%csfs_Ji(icsf_LS)%Ji(2)%      &
                                                              items(i2)
!write(iwrite_lsjj,*)'          J2  ',J2
                   J12_min=abs(Ji(1)-Ji(2))
                   J12_max=Ji(1)+Ji(2)
                   do J12=J12_min,J12_max,2
                     J_i(2)=J12
!write(iwrite_lsjj,*)'             J12 ',J12
                     if(J12.eq.iJ_total) then
                        icsf_JJ=icsf_JJ+1
                        if(icsf_JJ.gt.expansion_cfg_JJ%size) then
                           write(*,*)      &
                       'STOP at subroutine ... module transform_lsjj:',&
                          ' icsf_JJ.gt.expansion_cfg_JJ%size'
                           stop
                        end if
!allocate arrays of csf
                       expansion_cfg_JJ%csfs(icsf_JJ)%nosubc =         &
                                             cfg_Ji_structure%nr_of_subc
                       allocate(expansion_cfg_JJ%csfs(icsf_JJ)%        &
                        subc_cfg(expansion_cfg_JJ%csfs(icsf_JJ)%nosubc))
                       allocate(expansion_cfg_JJ%csfs(icsf_JJ)%        &
                        subc(expansion_cfg_JJ%csfs(icsf_JJ)%nosubc))
                       allocate(expansion_cfg_JJ%csfs(icsf_JJ)%iM1(    &
                        expansion_cfg_JJ%csfs(icsf_JJ)%nosubc))
                       allocate(expansion_cfg_JJ%csfs(icsf_JJ)%        &
                        iM2(expansion_cfg_JJ%csfs(icsf_JJ)%nosubc))
                       allocate(expansion_cfg_JJ%csfs(icsf_JJ)%        &
                        iJ(expansion_cfg_JJ%csfs(icsf_JJ)%nosubc))
!assign values
                       do isubc=1,expansion_cfg_JJ%csfs(icsf_JJ)%nosubc,1
                        expansion_cfg_JJ%csfs(icsf_JJ)%iJ(isubc) = 0
                        expansion_cfg_JJ%csfs(icsf_JJ)%subc_cfg(isubc) &
                         =expansion_cfg_LS%csfs(inr_of_csf_in_expansion_LS)%subc_cfg(isubc)
                        expansion_cfg_JJ%csfs(icsf_JJ)%subc(isubc) =  &
                         expansion_cfg_LS%csfs(inr_of_csf_in_expansion_LS)%subc(isubc)
                         expansion_cfg_JJ%csfs(icsf_JJ)%iM1(isubc)=Ji(isubc)
                         expansion_cfg_JJ%csfs(icsf_JJ)%iM2(isubc)=J_i(isubc)
                       end do
                     end if
                   end do
                 end do
               end if
            end do
         end do
!
         if(associated(Ji)) deallocate(Ji, STAT = error)
         if(associated(J_i)) deallocate(J_i, STAT = error)
!
      end if
!dealloacte cfg_Ji_structure
      do icsf=1,cfg_Ji_structure%nr_of_nomach_csfs
        do isubc=1,cfg_Ji_structure%nr_of_subc
         if(associated(cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%items))&
          deallocate(cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%items,STAT=error)
        end do
        if(associated(cfg_Ji_structure%csfs_Ji)) deallocate(cfg_Ji_structure%csfs_Ji(icsf)%Ji, STAT = error)
      end do
      if (associated(cfg_Ji_structure%csfs_Ji)) deallocate(cfg_Ji_structure%csfs_Ji, STAT = error)
!
!write(*,*) '      end subroutine form_csfs_JJ'
!!write(iwrite_log,*) '      end subroutine form_csfs_JJ'
!
      contains
!
!--- 	 subroutine  form_Ji_structure  ---------
         subroutine form_Ji_structure
         implicit none
!items in list shows the serial numbers of nonmaching
!csfs_LS in state%list_of_csfs_LS
         integer::iJ
!integer::inr_of_csf_in_state
!integer::inr_of_csf_in_csfs_LS
         integer::icsf, isubc
         integer::iSi, iLi, iJi_min, iJi_max
         integer::inumber_of_Ji, item, inr_of_subc
         integer::icsf_nr_in_expansion_LS
         inr_of_subc=expansion_cfg_LS%csfs(1)%nosubc
!write(iwrite_cfg_expansions_JJ,*)'   nr of nonmaching csf: ',nomach_csfs_LS%list_size
         cfg_Ji_structure%nr_of_subc = inr_of_subc
         cfg_Ji_structure%nr_of_nomach_csfs = nomach_csfs_LS%list_size
!cfg_Ji_structure%nr_of_nomach_csfs = inomach
         allocate(cfg_Ji_structure%csfs_Ji(cfg_Ji_structure%nr_of_nomach_csfs))
!toliau nustatom ir alokuojam csfs_Ji turini
!1. for each nomach_csf allocate	Ji list
         do icsf=1,cfg_Ji_structure%nr_of_nomach_csfs
            allocate(cfg_Ji_structure%csfs_Ji(icsf)%Ji(cfg_Ji_structure%nr_of_subc))
!cfg_Ji_structure%csfs_Ji(icsf)%nr_of_csf=inr_of_csf_in_csfs_LS
            icsf_nr_in_expansion_LS=nomach_csfs_LS%items(icsf)
            cfg_Ji_structure%csfs_Ji(icsf)%nr_of_csf=icsf_nr_in_expansion_LS
!2. fill in Ji list
            do isubc=1,cfg_Ji_structure%nr_of_subc
               iLi=expansion_cfg_LS%csfs(icsf_nr_in_expansion_LS)%subc(isubc)%iL
               iSi=expansion_cfg_LS%csfs(icsf_nr_in_expansion_LS)%subc(isubc)%iS
!count J_min, J_max for each subshell
!define possible Ji values and store them
               iJi_min=abs(iLi-iSi)
               iJi_max=abs(iLi+iSi)
               inumber_of_Ji = int((iJi_max-iJi_min)/2)+1
               cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%list_size=inumber_of_Ji
               allocate(cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%items(cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%list_size))
               iJ=iJi_min
               do item=1,cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%list_size
                  if(iJ.gt.iJi_max) stop 'STOP at subroutine ... module transform_lsll: iJ.gt.iJi_max'
                  cfg_Ji_structure%csfs_Ji(icsf)%Ji(isubc)%items(item) = iJ
                  iJ=iJ+2
               end do
            end do
         end do
!deallocate nomach_csfs_LS%items
!they was needed for Ji_structure formation
         end subroutine form_Ji_structure
!
!--- 	 subroutine  define_number_of_csfs_JJ  ---------
!
         subroutine  define_number_of_csfs_JJ(irez)
!--------------------------------------------------------------------
! This subroutine defines the number of csfs in
! JJ coupling
!--------------------------------------------------------------------
         implicit none
         integer, intent(out)::irez
         integer::J1, J2, J12_min,J12_max,J12
         integer::isubc_aviable
         integer::icsf, i1, i2
         isubc_aviable = 2
         if(cfg_Ji_structure%nr_of_subc.gt.isubc_aviable) then
            write(*,*) "Number of subshells is available =", &
            isubc_aviable, &
            " we have =",expansion_cfg_LS%csfs(1)%nosubc
            write(*,*) "STOP at subroutine define_number_of_csfs_JJ ", &
            "module transform_lsjj: ",&
            "cfg_Ji_structure%nr_of_subc.gt.isubc_aviable"
            stop
         end if
         irez=0
         do icsf=1,cfg_Ji_structure%nr_of_nomach_csfs
            do i1=1,cfg_Ji_structure%csfs_Ji(icsf)%Ji(1)%list_size
               J1=cfg_Ji_structure%csfs_Ji(icsf)%Ji(1)%items(i1)
               if(cfg_Ji_structure%nr_of_subc.gt.1) then
                  do i2=1,cfg_Ji_structure%csfs_Ji(icsf)%Ji(2)%list_size
                     J2=cfg_Ji_structure%csfs_Ji(icsf)%Ji(2)%items(i2)
                     J12_min=abs(J1-J2)
                     J12_max=J1+J2
                     do J12=J12_min,J12_max,2
                        if(J12.eq.iJ_total) then
                           irez=irez+1
                        end if
                     end do
                  end do
               end if
            end do
         end do
         end subroutine  define_number_of_csfs_JJ
      end subroutine form_csfs_JJ
!
!--- 	 subroutine  matrix_LS_JJ  ---------
!
      subroutine matrix_LS_JJ(icsf_LS,icsf_JJ,rez)
!--------------------------------------------------------------------
! This subroutine calculates the transformation
! matrix element between two csfs
! (in LS and JJ couplings)
!--------------------------------------------------------------------
      implicit none
      integer,intent(in):: icsf_LS,icsf_JJ
      real(kind=dp),intent(out)::rez
      real(kind=dp)::rez_2
!real(kind=dp)::ninej, multp
      integer::isubc
      integer::iLi,iSi,iS_i,iL_i,iL_im1,iS_im1,iJi,iJ_i,iJ_im1
      integer(kind=i4b)::i,inn
      integer::nosubc !number of subshells
!
      i=2
      inn=2
!
      rez=ZERO_dp
      if(the_same_subc(expansion_cfg_LS%csfs(icsf_LS),expansion_cfg_JJ%csfs(icsf_JJ)))then
         nosubc=expansion_cfg_LS%csfs(icsf_LS)%nosubc
         rez=ONE_dp
         if(nosubc.gt.1) then
            do isubc=2,nosubc
               iLi=expansion_cfg_LS%csfs(icsf_LS)%subc(isubc)%iL
               iSi=expansion_cfg_LS%csfs(icsf_LS)%subc(isubc)%iS
               iL_i=expansion_cfg_LS%csfs(icsf_LS)%iM1(isubc)
               iS_i=expansion_cfg_LS%csfs(icsf_LS)%iM2(isubc)
               iL_im1=expansion_cfg_LS%csfs(icsf_LS)%iM1(isubc-1)
               iS_im1=expansion_cfg_LS%csfs(icsf_LS)%iM2(isubc-1)
               iJi=expansion_cfg_JJ%csfs(icsf_JJ)%iM1(isubc)
               iJ_i=expansion_cfg_JJ%csfs(icsf_JJ)%iM2(isubc)
               iJ_im1=expansion_cfg_JJ%csfs(icsf_JJ)%iM2(isubc-1)
!
               call JJPER(iL_im1,iLi,iS_im1,iSi,iL_i,iS_i,iJ_i,iJ_im1,iJi,rez_2)
!
!call NINE(iL_im1,iLi,iL_i,  &
!			 iS_im1,iSi,iS_i,  &
!			 iJ_im1,iJi,iJ_i,  &
!			 i,inn,ninej)
!			 ! last line for techn. and output parameters
!multp=dsqrt(dble((iL_i+1)*(iS_i+1)*(iJ_im1+1)*(iJi+1)))
!
               rez=rez*rez_2
!
            end do
         end if
      end if
!
      end subroutine matrix_LS_JJ
!
!
!--- 	 subroutine  dotransform_lsjj  ---------
!
      subroutine dotransform_lsjj
!--------------------------------------------------------------------
! This subroutine calculates the weights of the
! expansions in the JJ coupling
!--------------------------------------------------------------------
      implicit none
      real(kind=dp)::coeff_JJ, coeff_LS, melement, sum_of_state
      integer::icsf_LS,icsf_JJ
      sum_of_state = ZERO_dp
      do icsf_JJ=1,expansion_cfg_JJ%size
         coeff_JJ = ZERO_dp
         do icsf_LS=1, expansion_cfg_LS%size
            call matrix_LS_JJ(icsf_LS,icsf_JJ,melement)
            coeff_LS=expansion_cfg_LS%coeffs(icsf_LS)
            coeff_JJ = coeff_JJ + coeff_LS*melement
         end do !icsf_LS
!check whether |coeff_JJ| <= 1
         if((dabs(coeff_JJ)-dp_coeff_precision).gt.ONE_dp) then
            write(*,*)'possible error at subroutine dotransform_lsjj: coeff_JJ=',coeff_JJ
            write(iwrite_log,*)'possible error at subroutine dotransform_lsjj: coeff_JJ=',coeff_JJ
         end if
         expansion_cfg_JJ%coeffs(icsf_JJ)= coeff_JJ
         sum_of_state = sum_of_state + coeff_JJ*coeff_JJ
      end do ! icsf_JJ
!check whether |sum_of_state| <= 1
      if((sum_of_state-TWO_dp*dp_coeff_precision).gt.ONE_dp) then
         write(*,*)'possible error at subroutine dotransform_lsjj: sum_of_state=',sum_of_state
         write(iwrite_log,*)'possible error at subroutine dotransform_lsjj: sum_of_state=',sum_of_state
      end if
!write(iwrite_cfg_expansions_JJ,'(25x,f10.7)')sum_of_state
      end subroutine dotransform_lsjj
!
!------   subroutine print_cfg_lsjj      -------------------
!
      subroutine print_cfg_lsjj(itype)
!--------------------------------------------------------------------
! This subroutine prints the oneconfigurational
! expansions to the unit "iwrite_cfg_expansions_JJ"
! (see module "Coupling_constants")
!--------------------------------------------------------------------
      implicit none
      integer, intent(in) :: itype
      integer             :: icsf, isubc, j
      character(len=1) :: CVAL
      character(len=4) :: JVAL
!
!     LS - Coupling
!
      if(itype.gt.0) then
         write(iwrite_cfg_expansions_JJ,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JJ,*)                           &
         'state Nr.',expansion_cfg_LS%nr_of_state
         write(iwrite_cfg_expansions_JJ,'(3x,a3,a4,a9,3x,f15.8)')    &
         'J =', JVAL(states%states(expansion_cfg_LS%nr_of_state)%J), &
         ' Energy =',states%states(expansion_cfg_LS%nr_of_state)%energy
         write(iwrite_cfg_expansions_JJ,'(3x,a29,i2)')               &
         'expansion size (LS coupling): ', expansion_cfg_LS%size
         write(iwrite_cfg_expansions_JJ,*)''
         if(associated(expansion_cfg_LS%csfs)) then
            do icsf=1,expansion_cfg_LS%size
               if(associated(expansion_cfg_LS%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_JJ,*)'csf Nr.'
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JJ,                  &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf, &
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%in,      &
                   CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(1)%iN_big,  &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JJ,                  &
                    '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%in,     &
                   CVAL(1,expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_LS%csfs(icsf)%subc_cfg(j)%iN_big,  &
                     j=1,2),                                          &
                     'coeff:',expansion_cfg_LS%coeffs(icsf)
                  else
                    STOP &
                  'To many couled shells in Coupling_transform_cfg_LSJJ'
                  end if
               else
                  write(iwrite_cfg_expansions_JJ,'(9x,a51)')         &
                  'expansion_cfg_LS%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_LS%csfs(icsf)%subc)       &
                  .and.associated(expansion_cfg_LS%csfs(icsf)%iM1)   &
                  .and.associated(expansion_cfg_LS%csfs(icsf)%iM2)) then
                  if(expansion_cfg_LS%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JJ,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_LS%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_LS%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_LS%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JJ,                 &
                     '(12x,7(1x,i2,a1,i1))')                         &
                     (expansion_cfg_LS%csfs(icsf)%subc(j)%iS+1,      &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%subc(j)%iL), &
                     expansion_cfg_LS%csfs(icsf)%subc(j)%inr,        &
                     j=1,2),                                         &
                     expansion_cfg_LS%csfs(icsf)%iM2(2)+1,           &
                     CVAL(2,expansion_cfg_LS%csfs(icsf)%iM1(2))
                  end if
               else
                  write(iwrite_cfg_expansions_JJ,'(9x,a63)')&
       'expansion_cfg_LS%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JJ,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JJ,*)&
            'expansion_cfg_LS%csfs NOT associated'
         end if
      end if
!
!     JJ - Coupling
!
      if(itype.gt.1) then
         write(iwrite_cfg_expansions_JJ,*) &
         '-------------------------------------'
         write(iwrite_cfg_expansions_JJ,*) &
         'state Nr.',expansion_cfg_JJ%nr_of_state
         write(iwrite_cfg_expansions_JJ,'(3x,a3,a4,a9,3x,f15.8)')      &
        'J =', JVAL(states%states(expansion_cfg_JJ%nr_of_state)%J),    &
         ' Energy =',states%states(expansion_cfg_JJ%nr_of_state)%energy
         write(iwrite_cfg_expansions_JJ,'(3x,a29,i2)')                 &
         'expansion size (JJ coupling): ', expansion_cfg_JJ%size
         write(iwrite_cfg_expansions_JJ,*)''
         if(associated(expansion_cfg_JJ%csfs)) then
            do icsf=1,expansion_cfg_JJ%size
               if(associated(expansion_cfg_JJ%csfs(icsf)%subc_cfg)) then
                  if(icsf ==1)write(iwrite_cfg_expansions_JJ,*)'csf Nr.'
                  if(expansion_cfg_JJ%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JJ,                  &
                     '(1x,i5,3x,i2,a1,"(",i1,")",6x,a6,f10.7)') icsf, &
                      expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%in,     &
                   CVAL(1,expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%il),&
                     expansion_cfg_JJ%csfs(icsf)%subc_cfg(1)%iN_big,  &
                     'coeff:',expansion_cfg_JJ%coeffs(icsf)
                  else if (expansion_cfg_JJ%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JJ,                   &
                     '(1x,i5,3x,2(i2,a1,"(",i1,")"),6x,a6,f10.7)')icsf,&
                     (expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%in,      &
                    CVAL(1,expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%il),&
                     expansion_cfg_JJ%csfs(icsf)%subc_cfg(j)%iN_big,   &
                     j=1,2),                                           &
                     'coeff:',expansion_cfg_JJ%coeffs(icsf)
                  else
                    STOP &
                  'To many couled shells in Coupling_transform_cfg_LSJJ'
                  end if
               else
                  write(iwrite_cfg_expansions_JJ,'(9x,a51)') &
                  'expansion_cfg_JJ%csfs(icsf)%subc_cfg NOT associated'
               end if
               if(associated(expansion_cfg_JJ%csfs(icsf)%subc)       &
                  .and.associated(expansion_cfg_JJ%csfs(icsf)%iM1)   &
                  .and.associated(expansion_cfg_JJ%csfs(icsf)%iM2)) then
                  if(expansion_cfg_JJ%csfs(icsf)%nosubc == 1) then
                     write(iwrite_cfg_expansions_JJ,                &
                     '(13x,i2,a1,i1)')                              &
                     expansion_cfg_JJ%csfs(icsf)%subc(1)%iS+1,      &
                     CVAL(2,expansion_cfg_JJ%csfs(icsf)%subc(1)%iL),&
                     expansion_cfg_JJ%csfs(icsf)%subc(1)%inr
                  else if (expansion_cfg_JJ%csfs(icsf)%nosubc == 2) then
                     write(iwrite_cfg_expansions_JJ,                 &
                     '(12x,2(1x,i2,a1,i1),3x,"(",a4,","a4,")",a4)')  &
                     (expansion_cfg_JJ%csfs(icsf)%subc(j)%iS+1,      &
                     CVAL(2,expansion_cfg_JJ%csfs(icsf)%subc(j)%iL), &
                     expansion_cfg_JJ%csfs(icsf)%subc(j)%inr,        &
                     j=1,2),                                         &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM2(1)),       &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM1(2)),       &
                     JVAL(expansion_cfg_JJ%csfs(icsf)%iM2(2))
                  end if
               else
                  write(iwrite_cfg_expansions_JJ,'(9x,a63)') &
       'expansion_cfg_JJ%csfs(icsf)%subc (or iM1 or iM2) NOT associated'
               end if
               write(iwrite_cfg_expansions_JJ,*)' '
            end do
         else
            write(iwrite_cfg_expansions_JJ,*) &
            'expansion_cfg_JJ%csfs NOT associated'
         end if
         write(iwrite_cfg_expansions_JJ,*) &
         '-------------------------------------'
      end if
      end subroutine print_cfg_lsjj
      end module Coupling_transform_cfg_LSJJ
