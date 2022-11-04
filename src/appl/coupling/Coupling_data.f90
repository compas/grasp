!
!***********************************************************************
!                                                                      *
      module Coupling_data
!                                                                      *
!     Written by G. Gaigalas,                                          *
!                                                                      *
!                                                                      *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use Coupling_constants
      use Coupling_structures
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      public :: fill_up_coupling_descriptions
      integer, parameter, public :: NR_OF_AVIALABLE_COUPLINGS  = 13
!     coupling_descriptions
      type(coupling_description), dimension(NR_OF_AVIALABLE_COUPLINGS) &
                                                :: coupling_descriptions
!     all the states under consideration
      type(set_of_states)::states
!     all the expansions under consideration
      type(set_of_set_of_coupling_expansions) ::all_expansions
!     all the classification data under consideration
      type(couplings_classification_data) :: all_classifications
!     all the classifications evaluation data under consideration
      type(couplings_evaluation_data) :: all_evaluations
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine fill_up_coupling_descriptions
!                                                                      *
!     This subroutine fills the array coupling_descriptions            *
!     (the array is used for printing - short and long names           *
!     of the couplings)                                                *
!                                                                      *
!     1 - LS coupling                                                  *
!     2 - JJ coupling                                                  *
!     3 - LK coupling                                                  *
!     4 - JK coupling                                                  *
!     5 - LS3 coupling                                                 *
!     6 - LSJ3 coupling                                                *
!     7 - LK3 coupling                                                 *
!     8 - JK3 coupling                                                 *
!     9 - cLSJ3 coupling                                               *
!    10 - LScjj coupling                                               *
!    11 - jj1 coupling                                                 *
!    12 - jj2 coupling                                                 *
!    13 - jj3 coupling                                                 *
!                                                                      *
!***********************************************************************
!
      coupling_descriptions(1)%short_name='LS'
      coupling_descriptions(1)%long_name ='LS coupling'
      coupling_descriptions(1)%iM1_name  ='L_i'
      coupling_descriptions(1)%iM2_name  ='S_i'
!
      coupling_descriptions(2)%short_name='JJ'
      coupling_descriptions(2)%long_name ='JJ coupling'
      coupling_descriptions(2)%iM1_name  ='Ji'
      coupling_descriptions(2)%iM2_name  ='J_i'
!
      coupling_descriptions(3)%short_name='LK'
      coupling_descriptions(3)%long_name ='LK coupling'
      coupling_descriptions(3)%iM1_name  ='iM1'
      coupling_descriptions(3)%iM2_name  ='iM2'
!
      coupling_descriptions(4)%short_name='JK'
      coupling_descriptions(4)%long_name ='JK coupling'
      coupling_descriptions(4)%iM1_name  ='iM1'
      coupling_descriptions(4)%iM2_name  ='iM2'
!
      coupling_descriptions(5)%short_name='LS3'
      coupling_descriptions(5)%long_name ='LS3 coupling'
      coupling_descriptions(5)%iM1_name  ='L_i'
      coupling_descriptions(5)%iM2_name  ='S_i'
!
      coupling_descriptions(6)%short_name='LSJ3'
      coupling_descriptions(6)%long_name ='LSJ3 coupling'
      coupling_descriptions(6)%iM1_name  ='L_i'
      coupling_descriptions(6)%iM2_name  ='S_i'
      coupling_descriptions(6)%iJ_name   ='J_i'
!
      coupling_descriptions(7)%short_name='LK3'
      coupling_descriptions(7)%long_name ='LK3 coupling'
      coupling_descriptions(7)%iM1_name  ='iM1'
      coupling_descriptions(7)%iM2_name  ='iM2'
!
      coupling_descriptions(8)%short_name='JK3'
      coupling_descriptions(8)%long_name ='JK3 coupling'
      coupling_descriptions(8)%iM1_name  ='iM1'
      coupling_descriptions(8)%iM2_name  ='iM2'
      coupling_descriptions(8)%iJ_name   ='J_i'
!
      coupling_descriptions(9)%short_name='cLSJ3'
      coupling_descriptions(9)%long_name ='cLSJ3 coupling'
      coupling_descriptions(9)%iM1_name  ='L_i'
      coupling_descriptions(9)%iM2_name  ='S_i'
      coupling_descriptions(9)%iJ_name   ='J_i'
!
      coupling_descriptions(10)%short_name='LScjj'
      coupling_descriptions(10)%long_name ='LScjj coupling'
      coupling_descriptions(10)%iM1_name  ='L_i'
      coupling_descriptions(10)%iM2_name  ='S_i'
      coupling_descriptions(10)%iJ_name   ='J_i'
!
      coupling_descriptions(11)%short_name='jj1'
      coupling_descriptions(11)%long_name ='jj1 coupling'
      coupling_descriptions(11)%iM1_name  ='L_i'
      coupling_descriptions(11)%iM2_name  ='S_i'
      coupling_descriptions(11)%iJ_name   ='J_i'
!
      coupling_descriptions(12)%short_name='jj2'
      coupling_descriptions(12)%long_name ='jj2 coupling'
      coupling_descriptions(12)%iM1_name  ='L_1'
      coupling_descriptions(12)%iM2_name  ='S_1'
      coupling_descriptions(12)%iJ_name   ='J_i'
!
      coupling_descriptions(13)%short_name='jj3'
      coupling_descriptions(13)%long_name ='jj3 coupling'
      coupling_descriptions(13)%iM1_name  ='L_i'
      coupling_descriptions(13)%iM2_name  ='S_i'
      coupling_descriptions(13)%iJ_name   ='J_i'
      end subroutine fill_up_coupling_descriptions
!
      end module Coupling_data
