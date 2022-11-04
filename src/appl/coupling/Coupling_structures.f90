!***********************************************************************
!
      module Coupling_structures
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use Coupling_constants
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      private :: is_equal_subc_LS_cfg
      private :: is_equal_subc_LS_state
      private :: is_equal_csf_LS
      private :: set_subc_LS_cfg_equal
      private :: set_subc_LS_state_equal
      private :: set_csf_LS_equal        !for "empty" b only !!!
      public :: the_same_subc
!-----------------------------------------------
!   D e f i n i t i o n s
!-----------------------------------------------
      type :: subc_LS_cfg ! subshell's configuration
         integer :: in
         integer :: il
         integer :: iN_big
      end type subc_LS_cfg
!
      type :: subc_LS_state ! type - subshell's state
         integer :: inr
         integer :: iL
         integer :: iS
      end type subc_LS_state
!
      type :: coupling_description
!         character(len=2) :: short_name
!         character(len=20) :: long_name
         character(len=3) :: short_name
         character(len=20) :: long_name
         character(len=3) :: iM1_name
         character(len=3) :: iM2_name
         character(len=3) :: iJ_name
      end type coupling_description
!
!masyviniai tipai
!
      type :: csf_LS
         integer :: nosubc
         type(subc_LS_cfg), dimension(:),pointer:: subc_cfg
         type(subc_LS_state), dimension(:),pointer:: subc
         integer, dimension(:),pointer :: iM1
         integer, dimension(:),pointer :: iM2
         integer, dimension(:),pointer :: iJ
      end type csf_LS
!
      type::expansion
!serial number of state in states%states() or ...
         integer::nr_of_state
         integer::size
         real(kind=dp),dimension(:), pointer :: coeffs
         type(csf_LS), dimension(:), pointer :: csfs
         character(len=20) :: note20
      end type expansion
!
!state quantum numbers
      type::state
         real(kind=dp)::energy
         integer::J   ! 2*J value
         byte::parity
      end type state
!
!------------------------------------------------------------
!	sets
!
      type :: set_of_csfs_LS
         integer :: nr_of_csfs
         type(csf_LS), dimension(:), pointer :: csfs
      end type set_of_csfs_LS
!
      type:: set_of_states
         integer::nr_of_states  !number of states under consideration
         type(state),dimension(:), pointer :: states
      end type set_of_states
!
!expansions
!
      type:: set_of_expansions
         integer::nr_of_expansions  !number of expansions
         type(expansion),dimension(:), pointer :: expansions
      end type set_of_expansions
!
      type:: set_of_coupling_expansions
         integer::icoupling  ! coupling serial nr in coupling descriptions
         integer::nr_of_expansions  !number of expansions
         type(expansion),dimension(:), pointer :: expansions
      end type set_of_coupling_expansions
!
      type:: set_of_set_of_coupling_expansions
         integer::nr_of_coupling_expansions  !number of coupling expansions
         type(set_of_coupling_expansions),dimension(:), pointer ::     &
                                                     coupling_expansions
      end type set_of_set_of_coupling_expansions
!
!--------------------------------------------------------------------
!	lists	- additional types used for technical purposes
      type:: list
         integer::list_size  !number items in a list
         integer, dimension(:),pointer:: items !serial numbers of lists items
      end type list
!
      type:: numbered_list
         integer:: number
         type(list)::list
      end type numbered_list
!
      type :: list_of_coeffs
         integer::list_size
         real(kind=dp),dimension(:), pointer::coeffs
      end type list_of_coeffs
!
!------------  other   --------------------
      type :: cfg_structure
         integer::nr_of_cfgs
         type(list), dimension(:), pointer::csfs_lists
      end type cfg_structure
!
!------------  evaluation  and classification  ----------------
!
      type :: R_values
!integer:: icoupling ! coupling serial nr in coupling descriptions
         integer:: nr_of_J
         integer, dimension(:), pointer :: J
         real(kind=dp),dimension(:), pointer::RJ_values
         real(kind=dp) :: R_value
         real(kind=dp) :: P_value
      end type R_values
!
      type :: couplings_evaluation_data
         integer::nr_of_couplings
         type(R_values), dimension(:), pointer :: R
      end type couplings_evaluation_data
!
!classification
!
      type :: state_classification_data
         type(csf_LS):: csf
         real(kind=dp) :: coeff
         real(kind=dp) :: coeffs_dispersion
         integer:: coupling_nr
      end type state_classification_data
!
      type :: coupling_classification_data
!integer:: icoupling ! coupling serial nr in coupling descriptions
         type(state_classification_data), dimension(:), pointer ::     &
                                                                  states
         integer:: nr_of_states
         logical :: is_one_to_one
      end type coupling_classification_data
!
      type :: couplings_classification_data
         integer :: nr_of_couplings
         type(coupling_classification_data), dimension(:), pointer ::  &
                                                               couplings
      end type couplings_classification_data
!
! Define several derived data structures to keep information about subshell
! states and terms together.
!
      type, public :: subshell_term
         integer :: j             ! Angular momentum 2*j.
         integer :: Q             ! Subshell total quasispin 2*Q.
         integer :: nu            ! Seniority number.
         integer :: subshellJ     ! Subshell total angular momentum 2*J.
         integer :: Nr            ! State identifier Nr.
      end type subshell_term
!
      type, public  :: subshell_term_LS
         integer :: l_shell ! angular momentum l
         integer :: w       ! w
         integer :: Q       ! Subshell total quasispin 2*Q
         integer :: LL      ! Subshell total angular momentum 2*L
         integer :: S       ! Subshell total angular momentum 2*S
      end type subshell_term_LS
!
      type :: LS_jj_me
!         sequence
         integer       :: w, Q, L, S, J, Nm, Qm, Jm, Qp, Jp
         integer       :: factor
         integer(selected_int_kind(18)) :: nom, denom
      end type LS_jj_me
!
! Define the values of all LS-jj transformation coefficients
!
      integer, parameter ::                                            &
                       LS_jj_number_p3 =  10, LS_jj_number_p4  =    9, &
                       LS_jj_number_p5 =   2, LS_jj_number_p6  =    1, &
                       LS_jj_number_d3 =  65, LS_jj_number_d4  =  166, &
                       LS_jj_number_d5 = 184, LS_jj_number_d6  =  166, &
                       LS_jj_number_d7 =  65, LS_jj_number_d8  =   19, &
                       LS_jj_number_d9 =   2, LS_jj_number_d10 =    1, &
                       LS_jj_number_f3 = 216, LS_jj_number_f4  = 1210, &
                       LS_jj_number_f5 =3799, LS_jj_number_f6  = 7313, &
                       LS_jj_number_f7 =8003
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      character(len=4), public :: symbol_I_Numb
      integer, public :: I_Numb
!
!---------------------------------------------------------------
!
      interface operator(==)
         module procedure is_equal_subc_LS_cfg
         module procedure is_equal_subc_LS_state
!module procedure is_equal_cfg_LS
         module procedure is_equal_csf_LS
      end interface
!
      interface assignment(=)
         module procedure set_subc_LS_cfg_equal
         module procedure set_subc_LS_state_equal
         module procedure set_csf_LS_equal        !for "empty" b only !!!
!module procedure set_expansions_equal
      end interface
!
contains
!
!***********************************************************************
!                                                                      *
      function is_equal_subc_LS_cfg(a,b)                    result(yes)
!                                                                      *
!     This function "defines" the equivalence of two variables of      *
!     type(subc_LS_cfg).                                               *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(subc_LS_cfg), intent(in) :: a, b
      logical                       :: yes
      if (a%in==b%in.and.a%il==b%il.and.a%iN_big==b%iN_big) then
         yes = .true.
      else
         yes = .false.
      end if
      end function is_equal_subc_LS_cfg
!
!***********************************************************************
!                                                                      *
      function is_equal_subc_LS_state(a,b)                  result(yes)
!                                                                      *
!     This function "defines" the equivalence of two variables of      *
!     type(subc_LS_state).                                             *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(subc_LS_state), intent(in) :: a, b
      logical                         :: yes
      if (a%inr==b%inr.and.a%iL==b%iL.and.a%iS==b%iS) then
!if (a%cfg==b%cfg.and.a%inr==b%inr.and.a%iL==b%iL.and.a%iS==b%iS) then
         yes = .true.
      else
         yes = .false.
      end if
      end function is_equal_subc_LS_state
!
!***********************************************************************
!                                                                      *
      function is_equal_csf_LS(a,b)                         result(yes)
!                                                                      *
!     This function "defines" the equivalence of two variables of      *
!     type(csf_LS).                                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(csf_LS), intent(in) :: a, b
      logical                  :: yes
      integer                  :: isubc
      if(a%nosubc.eq.b%nosubc) then
         yes=.true.
         if(a%nosubc.gt.0) then
            do isubc=1, a%nosubc
               if(.not. a%subc_cfg(isubc)==b%subc_cfg(isubc)) then
                  yes = .false.
                  exit
               else if(.not. a%subc(isubc)==b%subc(isubc)) then
                  yes = .false.
                  exit
               else if(.not. a%iM1(isubc).eq.b%iM1(isubc)) then
                  yes = .false.
                  exit
               else if(.not. a%iM2(isubc).eq.b%iM2(isubc)) then
                  yes = .false.
                  exit
               end if
            end do
         end if
      end if
      end function is_equal_csf_LS
!
!***********************************************************************
!                                                                      *
      subroutine set_subc_LS_cfg_equal (b, a)
!                                                                      *
!     This function "overloads" = for two variables of                 *
!     type(subc_LS_cfg).                                               *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      type(subc_LS_cfg), intent(in)  :: a
      type(subc_LS_cfg), intent(out) :: b
      b%in = a%in
      b%il = a%il
      b%iN_big = a%iN_big
      end  subroutine set_subc_LS_cfg_equal
!
!***********************************************************************
!                                                                      *
      subroutine set_subc_LS_state_equal (b, a)
!                                                                      *
!     This function "overloads" = for two variables of                 *
!     type(subc_LS_cfg).                                               *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      type(subc_LS_state), intent(in)  :: a
      type(subc_LS_state), intent(out) :: b
      b%inr = a%inr
      b%iL = a%iL
      b%iS = a%iS
      end  subroutine set_subc_LS_state_equal
!
!***********************************************************************
!                                                                      *
      subroutine set_csf_LS_equal (b, a)
!                                                                      *
!     This function "overloads" = for two variables of                 *
!     type(csf_LS).                                                    *
!                                                                      *
!     IN THIS PROGRAM VERSION SUBROTINE set_csf_LS_equal (b, a)        *
!     AVIABLE ONLY FOR THE CASE OF "EMPTY" ARGUMENT "b"                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      type(csf_LS), intent(in)  :: a
      type(csf_LS), intent(out) :: b
      b%nosubc = a%nosubc
      if(associated(a%subc).and.associated(a%subc_cfg).and.            &
         associated(a%iM1) .and.associated(a%iM2) .and.                &
                                                 associated(a%iJ)) then
         allocate(b%subc(b%nosubc))
         allocate(b%subc_cfg(b%nosubc))
         allocate(b%iM1(b%nosubc))
         allocate(b%iM2(b%nosubc))
         allocate(b%iJ(b%nosubc))
         do i=1, b%nosubc
            b%subc(i)=a%subc(i)
            b%subc_cfg(i)=a%subc_cfg(i)
            b%iM1(i)=a%iM1(i)
            b%iM2(i)=a%iM2(i)
            b%iJ(i)=a%iJ(i)
         end do
      end if
      end  subroutine set_csf_LS_equal
!
!***********************************************************************
!                                                                      *
      function the_same_subc(csf1, csf2)                    result(rez)
!                                                                      *
!     This function checks whether two csfs have                       *
!     the same configuration and states of subshells                   *
!     (i.e. differ intermediate momentas only)                         *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(csf_LS)::csf1, csf2
      integer::isubc
      logical::rez
      rez=.true.
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            else if(.not.csf1%subc(isubc)==csf2%subc(isubc)) then
               rez=.false.
               exit
            end if
         end do
      else
         rez=.false.
      end if
      end function the_same_subc
!
!***********************************************************************
!                                                                      *
      function the_same_cfg(csf1, csf2)                     result(rez)
!                                                                      *
!     This function checks whether two csfs is of the same             *
!     configuration                                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(csf_LS)::csf1, csf2
      integer::isubc
      logical::rez
      rez=.true.
      if(csf1%nosubc.eq.csf2%nosubc) then
         do isubc=1,csf1%nosubc
            if(.not.csf1%subc_cfg(isubc)==csf2%subc_cfg(isubc)) then
               rez=.false.
               exit
            end if
         end do
      else
         rez=.false.
      end if
      end function the_same_cfg
!
!***********************************************************************
!                                                                      *
      subroutine print_list(iwrite, inputlist)
!                                                                      *
!     This subroutine prints list "inputlist" to the unit "iwrite"     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
      type(list), intent(in)::inputlist
      integer, intent(in) :: iwrite
      integer::i
      write(iwrite,*) 'subroutine print_list:'
      write(iwrite,*) 'list_size:', inputlist%list_size
      if(associated(inputlist%items).and.inputlist%list_size.gt.0) then
         do i=1, inputlist%list_size
            write(iwrite, *)i, inputlist%items(i)
         end do
      end if
      write(iwrite,*) 'end subroutine print_list:'
      end subroutine
!
      end module Coupling_structures
