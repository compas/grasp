!***********************************************************************
!                                                                      *
      module Coupling_inside_shell
!
!-----------------------------------------------------------------------
!                                                                      *
!    This module contains the procedures which are used to control     *
!    the jj-LS transformation insade the shell of equivalents          *
!    electrons                                                         *
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
      USE Coupling_structures
      USE jj2lsj_data_1_C
      USE jj2lsj_data_2_C
      USE jj2lsj_data_3_C
      implicit none
!
      public   :: gettermjj
                 ! This procedure return all allowed subshell terms
                 ! (j, w, Q, J) for given j^N which must be 1/2, 3/2,
                 ! 5/2, 7/2 or 9/2.
      public  :: gettermLS
!               This procedure return all allowed subshell terms
!               (l, w, Q, L, S) for given l^N which must be 0, 1, 2 or 3.
      public  :: gettermLSQ
!
      public   :: coefLSjj
                 ! Returns the value of the LS-jj transformation matrix
                 ! for a given set of quantum numbers.
      private  :: coefLSjj2
                 ! Returns the value of the LS-jj transformation matrix
                 ! (l^2 LSJ| j_1 j_2 J)
      public  :: coefLSjjs
                 ! Returns the value of the LS-jj transformation matrix
                 ! for the list.
      public  :: JTHN
      public  :: getJMinMax
      public  :: getjjJMinMax
!
! Define all possible antisymetric subshell terms in jj-coupling
!
   type(subshell_term), dimension(1:63), parameter, public ::      &
      terms_jj = (/                                                &
      !
      ! j = 1/2; these terms have indices terms_jj(1:2)
      subshell_term(1, 0, 1, 1, 0), subshell_term(1, 1, 0, 0, 0),  &
      !
      ! j = 3/2; these terms have indices terms_jj(3:5)
      subshell_term(3, 1, 1, 3, 0), subshell_term(3, 2, 0, 0, 0),  &
      subshell_term(3, 0, 2, 4, 0),                                &
      !
      ! j = 5/2; these terms have indices terms_jj(6:11)
      subshell_term(5, 2, 1, 5, 0), subshell_term(5, 0, 3, 3, 0),  &
      subshell_term(5, 0, 3, 9, 0), subshell_term(5, 3, 0, 0, 0),  &
      subshell_term(5, 1, 2, 4, 0), subshell_term(5, 1, 2, 8, 0),  &
      !
      ! j = 7/2; these terms have indices terms_jj(12:25)
      subshell_term(7, 1, 3, 3, 0), subshell_term(7, 1, 3, 5, 0),  &
      subshell_term(7, 3, 1, 7, 0), subshell_term(7, 1, 3, 9, 0),  &
      subshell_term(7, 1, 3,11, 0), subshell_term(7, 1, 3,15, 0),  &
      subshell_term(7, 4, 0, 0, 0), subshell_term(7, 2, 2, 4, 0),  &
      subshell_term(7, 0, 4, 4, 0), subshell_term(7, 2, 2, 8, 0),  &
      subshell_term(7, 0, 4, 8, 0), subshell_term(7, 0, 4,10, 0),  &
      subshell_term(7, 2, 2,12, 0), subshell_term(7, 0, 4,16, 0),  &
      !
      ! j = 9/2; these terms have indices terms_jj(26:63)
      subshell_term(9, 0, 5, 1, 0), subshell_term(9, 2, 3, 3, 0),  &
      subshell_term(9, 2, 3, 5, 0), subshell_term(9, 0, 5, 5, 0),  &
      subshell_term(9, 2, 3, 7, 0), subshell_term(9, 0, 5, 7, 0),  &
      subshell_term(9, 4, 1, 9, 0), subshell_term(9, 2, 3, 9, 0),  &
      subshell_term(9, 0, 5, 9, 0), subshell_term(9, 2, 3,11, 0),  &
      subshell_term(9, 0, 5,11, 0), subshell_term(9, 2, 3,13, 0),  &
      subshell_term(9, 0, 5,13, 0), subshell_term(9, 2, 3,15, 0),  &
      subshell_term(9, 0, 5,15, 0), subshell_term(9, 2, 3,17, 0),  &
      subshell_term(9, 0, 5,17, 0), subshell_term(9, 0, 5,19, 0),  &
      subshell_term(9, 2, 3,21, 0), subshell_term(9, 0, 5,25, 0),  &
      subshell_term(9, 5, 0, 0, 0), subshell_term(9, 1, 4, 0, 0),  &
      subshell_term(9, 3, 2, 4, 0), subshell_term(9, 1, 4, 4, 0),  &
      subshell_term(9, 1, 4, 6, 0), subshell_term(9, 3, 2, 8, 0),  &
      subshell_term(9, 1, 4, 8, 1), subshell_term(9, 1, 4, 8, 2),  &
      subshell_term(9, 1, 4,10, 0), subshell_term(9, 3, 2,12, 0),  &
      subshell_term(9, 1, 4,12, 1), subshell_term(9, 1, 4,12, 2),  &
      subshell_term(9, 1, 4,14, 0), subshell_term(9, 3, 2,16, 0),  &
      subshell_term(9, 1, 4,16, 0), subshell_term(9, 1, 4,18, 0),  &
      subshell_term(9, 1, 4,20, 0), subshell_term(9, 1, 4,24, 0) /)
      !
      ! j = 11/2;
   type(subshell_term), dimension(1:1), parameter ::   &
      terms_jj_11_1 =(/                                &
      subshell_term(11, 5, 1,11, 0)  /)
   type(subshell_term), dimension(1:6), parameter ::   &
      terms_jj_11_2 =(/                                &
      subshell_term(11, 6, 0, 0, 0),                   &
      subshell_term(11, 4, 2, 4, 0),                   &
      subshell_term(11, 4, 2, 8, 0),                   &
      subshell_term(11, 4, 2,12, 0),                   &
      subshell_term(11, 4, 2,16, 0),                   &
      subshell_term(11, 4, 2,20, 0)  /)
      !
      ! j = 13/2;
   type(subshell_term), dimension(1:1), parameter ::   &
      terms_jj_13_1 =(/                                &
      subshell_term(13, 6, 1,13, 0)  /)
   type(subshell_term), dimension(1:7), parameter ::   &
      terms_jj_13_2 =(/                                &
      subshell_term(13, 7, 0, 0, 0),                   &
      subshell_term(13, 5, 2, 4, 0),                   &
      subshell_term(13, 5, 2, 8, 0),                   &
      subshell_term(13, 5, 2,12, 0),                   &
      subshell_term(13, 5, 2,16, 0),                   &
      subshell_term(13, 5, 2,20, 0),                   &
      subshell_term(13, 5, 2,24, 0)  /)
      !
      ! j = 15/2;
   type(subshell_term), dimension(1:1), parameter ::   &
      terms_jj_15_1 =(/                                &
      subshell_term(15, 7, 1,15, 0)  /)
   type(subshell_term), dimension(1:8), parameter ::   &
      terms_jj_15_2 =(/                                &
      subshell_term(15, 8, 0, 0, 0),                   &
      subshell_term(15, 6, 2, 4, 0),                   &
      subshell_term(15, 6, 2, 8, 0),                   &
      subshell_term(15, 6, 2,12, 0),                   &
      subshell_term(15, 6, 2,16, 0),                   &
      subshell_term(15, 6, 2,20, 0),                   &
      subshell_term(15, 6, 2,24, 0),                   &
      subshell_term(15, 6, 2,28, 0)  /)
   !
   ! Define all possible antisymetric subshell terms in LS-coupling
   !
   type(subshell_term_LS), dimension(1:2), parameter, public ::   &
      term_LS_s =(/                                               &
      subshell_term_LS(0, 0, 0, 0, 1),                            &
      subshell_term_LS(0, 0, 1, 0, 0)  /)
   type(subshell_term_LS), dimension(1:6), parameter, public ::   &
      term_LS_p =(/                                               &
      subshell_term_LS(1, 0, 0, 0, 3),                            &
      subshell_term_LS(1, 0, 2, 2, 1),                            &
      subshell_term_LS(1, 0, 0, 4, 1),                            &
      subshell_term_LS(1, 0, 3, 0, 0),                            &
      subshell_term_LS(1, 0, 1, 2, 2),                            &
      subshell_term_LS(1, 0, 1, 4, 0)  /)
   type(subshell_term_LS), dimension(1:32), parameter, public ::  &
      term_LS_d =(/                                               &
      subshell_term_LS(2, 0, 0, 0, 5),                            &
      subshell_term_LS(2, 0, 0, 0, 1),                            &
      subshell_term_LS(2, 0, 2, 2, 3),                            &
      subshell_term_LS(2, 0, 2, 2, 1),                            &
      subshell_term_LS(2, 0, 4, 4, 1),                            &
      subshell_term_LS(2, 0, 2, 4, 1),                            &
      subshell_term_LS(2, 0, 0, 4, 3),                            &
      subshell_term_LS(2, 0, 0, 4, 1),                            &
      subshell_term_LS(2, 0, 2, 6, 3),                            &
      subshell_term_LS(2, 0, 2, 6, 1),                            &
      subshell_term_LS(2, 0, 0, 6, 1),                            &
      subshell_term_LS(2, 0, 2, 8, 1),                            &
      subshell_term_LS(2, 0, 0, 8, 3),                            &
      subshell_term_LS(2, 0, 0, 8, 1),                            &
      subshell_term_LS(2, 0, 2,10, 1),                            &
      subshell_term_LS(2, 0, 0,12, 1),                            &
      subshell_term_LS(2, 0, 5, 0, 0),                            &
      subshell_term_LS(2, 0, 1, 0, 0),                            &
      subshell_term_LS(2, 0, 3, 2, 2),                            &
      subshell_term_LS(2, 0, 1, 2, 2),                            &
      subshell_term_LS(2, 0, 1, 4, 4),                            &
      subshell_term_LS(2, 0, 1, 4, 2),                            &
      subshell_term_LS(2, 0, 3, 4, 0),                            &
      subshell_term_LS(2, 0, 1, 4, 0),                            &
      subshell_term_LS(2, 0, 3, 6, 2),                            &
      subshell_term_LS(2, 0, 1, 6, 2),                            &
      subshell_term_LS(2, 0, 1, 6, 0),                            &
      subshell_term_LS(2, 0, 1, 8, 2),                            &
      subshell_term_LS(2, 0, 3, 8, 0),                            &
      subshell_term_LS(2, 0, 1, 8, 0),                            &
      subshell_term_LS(2, 0, 1,10, 2),                            &
      subshell_term_LS(2, 0, 1,12, 0)  /)
   type(subshell_term_LS), dimension(1:238), parameter, public :: &
      term_LS_f =(/                                               &
      subshell_term_LS(3, 0, 0, 0, 7),                            &
      subshell_term_LS(3, 0, 2, 2, 5),                            &
      subshell_term_LS(3, 0, 0, 4, 5),                            &
      subshell_term_LS(3, 0, 2, 6, 5),                            &
      subshell_term_LS(3, 0, 0, 8, 5),                            &
      subshell_term_LS(3, 0, 2,10, 5),                            &
      subshell_term_LS(3, 0, 0,12, 5),                            &
      subshell_term_LS(3, 1, 4, 0, 3),                            &
      subshell_term_LS(3, 2, 0, 0, 3),                            &
      subshell_term_LS(3, 1, 2, 2, 3),                            &
      subshell_term_LS(3, 2, 2, 2, 3),                            &
      subshell_term_LS(3, 1, 4, 4, 3),                            &
      subshell_term_LS(3, 2, 2, 4, 3),                            &
      subshell_term_LS(3, 3, 2, 4, 3),                            &
      subshell_term_LS(3, 4, 0, 4, 3),                            &
      subshell_term_LS(3, 5, 0, 4, 3),                            &
      subshell_term_LS(3, 6, 0, 4, 3),                            &
      subshell_term_LS(3, 1, 4, 6, 3),                            &
      subshell_term_LS(3, 2, 2, 6, 3),                            &
      subshell_term_LS(3, 3, 2, 6, 3),                            &
      subshell_term_LS(3, 4, 2, 6, 3),                            &
      subshell_term_LS(3, 5, 0, 6, 3),                            &
      subshell_term_LS(3, 1, 4, 8, 3),                            &
      subshell_term_LS(3, 2, 2, 8, 3),                            &
      subshell_term_LS(3, 3, 2, 8, 3),                            &
      subshell_term_LS(3, 4, 2, 8, 3),                            &
      subshell_term_LS(3, 5, 0, 8, 3),                            &
      subshell_term_LS(3, 6, 0, 8, 3),                            &
      subshell_term_LS(3, 7, 0, 8, 3),                            &
      subshell_term_LS(3, 1, 2,10, 3),                            &
      subshell_term_LS(3, 2, 2,10, 3),                            &
      subshell_term_LS(3, 3, 2,10, 3),                            &
      subshell_term_LS(3, 4, 0,10, 3),                            &
      subshell_term_LS(3, 5, 0,10, 3),                            &
      subshell_term_LS(3, 1, 4,12, 3),                            &
      subshell_term_LS(3, 2, 2,12, 3),                            &
      subshell_term_LS(3, 3, 2,12, 3),                            &
      subshell_term_LS(3, 4, 0,12, 3),                            &
      subshell_term_LS(3, 5, 0,12, 3),                            &
      subshell_term_LS(3, 1, 2,14, 3),                            &
      subshell_term_LS(3, 2, 2,14, 3),                            &
      subshell_term_LS(3, 3, 0,14, 3),                            &
      subshell_term_LS(3, 1, 2,16, 3),                            &
      subshell_term_LS(3, 2, 0,16, 3),                            &
      subshell_term_LS(3, 3, 0,16, 3),                            &
      subshell_term_LS(3, 0, 2,18, 3),                            &
      subshell_term_LS(3, 0, 0,20, 3),                            &
      subshell_term_LS(3, 1, 0, 0, 1),                            &
      subshell_term_LS(3, 2, 0, 0, 1),                            &
      subshell_term_LS(3, 1, 4, 2, 1),                            &
      subshell_term_LS(3, 2, 2, 2, 1),                            &
      subshell_term_LS(3, 3, 2, 2, 1),                            &
      subshell_term_LS(3, 4, 2, 2, 1),                            &
      subshell_term_LS(3, 5, 0, 2, 1),                            &
      subshell_term_LS(3, 1, 4, 4, 1),                            &
      subshell_term_LS(3, 2, 4, 4, 1),                            &
      subshell_term_LS(3, 3, 2, 4, 1),                            &
      subshell_term_LS(3, 4, 2, 4, 1),                            &
      subshell_term_LS(3, 5, 2, 4, 1),                            &
      subshell_term_LS(3, 6, 0, 4, 1),                            &
      subshell_term_LS(3, 7, 0, 4, 1),                            &
      subshell_term_LS(3, 1, 6, 6, 1),                            &
      subshell_term_LS(3, 2, 4, 6, 1),                            &
      subshell_term_LS(3, 3, 2, 6, 1),                            &
      subshell_term_LS(3, 4, 2, 6, 1),                            &
      subshell_term_LS(3, 5, 2, 6, 1),                            &
      subshell_term_LS(3, 6, 2, 6, 1),                            &
      subshell_term_LS(3, 7, 2, 6, 1),                            &
      subshell_term_LS(3, 8, 0, 6, 1),                            &
      subshell_term_LS(3, 9, 0, 6, 1),                            &
      subshell_term_LS(3,10, 0, 6, 1),                            &
      subshell_term_LS(3, 1, 4, 8, 1),                            &
      subshell_term_LS(3, 2, 4, 8, 1),                            &
      subshell_term_LS(3, 3, 2, 8, 1),                            &
      subshell_term_LS(3, 4, 2, 8, 1),                            &
      subshell_term_LS(3, 5, 2, 8, 1),                            &
      subshell_term_LS(3, 6, 2, 8, 1),                            &
      subshell_term_LS(3, 7, 0, 8, 1),                            &
      subshell_term_LS(3, 8, 0, 8, 1),                            &
      subshell_term_LS(3, 9, 0, 8, 1),                            &
      subshell_term_LS(3,10, 0, 8, 1),                            &
      subshell_term_LS(3, 1, 4,10, 1),                            &
      subshell_term_LS(3, 2, 4,10, 1),                            &
      subshell_term_LS(3, 3, 2,10, 1),                            &
      subshell_term_LS(3, 4, 2,10, 1),                            &
      subshell_term_LS(3, 5, 2,10, 1),                            &
      subshell_term_LS(3, 6, 2,10, 1),                            &
      subshell_term_LS(3, 7, 2,10, 1),                            &
      subshell_term_LS(3, 8, 0,10, 1),                            &
      subshell_term_LS(3, 9, 0,10, 1),                            &
      subshell_term_LS(3, 1, 4,12, 1),                            &
      subshell_term_LS(3, 2, 2,12, 1),                            &
      subshell_term_LS(3, 3, 2,12, 1),                            &
      subshell_term_LS(3, 4, 2,12, 1),                            &
      subshell_term_LS(3, 5, 2,12, 1),                            &
      subshell_term_LS(3, 6, 0,12, 1),                            &
      subshell_term_LS(3, 7, 0,12, 1),                            &
      subshell_term_LS(3, 8, 0,12, 1),                            &
      subshell_term_LS(3, 9, 0,12, 1),                            &
      subshell_term_LS(3, 1, 4,14, 1),                            &
      subshell_term_LS(3, 2, 2,14, 1),                            &
      subshell_term_LS(3, 3, 2,14, 1),                            &
      subshell_term_LS(3, 4, 2,14, 1),                            &
      subshell_term_LS(3, 5, 2,14, 1),                            &
      subshell_term_LS(3, 6, 0,14, 1),                            &
      subshell_term_LS(3, 7, 0,14, 1),                            &
      subshell_term_LS(3, 1, 4,16, 1),                            &
      subshell_term_LS(3, 2, 2,16, 1),                            &
      subshell_term_LS(3, 3, 2,16, 1),                            &
      subshell_term_LS(3, 4, 0,16, 1),                            &
      subshell_term_LS(3, 5, 0,16, 1),                            &
      subshell_term_LS(3, 1, 2,18, 1),                            &
      subshell_term_LS(3, 2, 2,18, 1),                            &
      subshell_term_LS(3, 3, 0,18, 1),                            &
      subshell_term_LS(3, 4, 0,18, 1),                            &
      subshell_term_LS(3, 1, 2,20, 1),                            &
      subshell_term_LS(3, 2, 0,20, 1),                            &
      subshell_term_LS(3, 0, 2,22, 1),                            &
      subshell_term_LS(3, 0, 0,24, 1),                            &
      subshell_term_LS(3, 0, 1, 6, 6),                            &
      subshell_term_LS(3, 1, 3, 4, 4),                            &
      subshell_term_LS(3, 2, 1, 4, 4),                            &
      subshell_term_LS(3, 3, 1, 4, 4),                            &
      subshell_term_LS(3, 1, 3, 6, 4),                            &
      subshell_term_LS(3, 2, 1, 6, 4),                            &
      subshell_term_LS(3, 1, 3, 8, 4),                            &
      subshell_term_LS(3, 2, 1, 8, 4),                            &
      subshell_term_LS(3, 3, 1, 8, 4),                            &
      subshell_term_LS(3, 0, 1, 2, 4),                            &
      subshell_term_LS(3, 1, 1,10, 4),                            &
      subshell_term_LS(3, 2, 1,10, 4),                            &
      subshell_term_LS(3, 0, 3, 0, 4),                            &
      subshell_term_LS(3, 1, 3,12, 4),                            &
      subshell_term_LS(3, 2, 1,12, 4),                            &
      subshell_term_LS(3, 0, 1,14, 4),                            &
      subshell_term_LS(3, 0, 1,16, 4),                            &
      subshell_term_LS(3, 1, 5, 6, 2),                            &
      subshell_term_LS(3, 2, 3, 6, 2),                            &
      subshell_term_LS(3, 6, 1, 6, 2),                            &
      subshell_term_LS(3, 8, 1, 6, 2),                            &
      subshell_term_LS(3, 1, 3, 4, 2),                            &
      subshell_term_LS(3, 2, 3, 4, 2),                            &
      subshell_term_LS(3, 3, 1, 4, 2),                            &
      subshell_term_LS(3, 4, 1, 4, 2),                            &
      subshell_term_LS(3, 3, 3, 6, 2),                            &
      subshell_term_LS(3, 5, 1, 6, 2),                            &
      subshell_term_LS(3, 1, 3, 8, 2),                            &
      subshell_term_LS(3, 2, 3, 8, 2),                            &
      subshell_term_LS(3, 4, 1, 8, 2),                            &
      subshell_term_LS(3, 5, 1, 8, 2),                            &
      subshell_term_LS(3, 5, 1, 4, 2),                            &
      subshell_term_LS(3, 4, 3, 6, 2),                            &
      subshell_term_LS(3, 7, 1, 6, 2),                            &
      subshell_term_LS(3, 9, 1, 6, 2),                            &
      subshell_term_LS(3, 3, 3, 8, 2),                            &
      subshell_term_LS(3, 6, 1, 8, 2),                            &
      subshell_term_LS(3, 7, 1, 8, 2),                            &
      subshell_term_LS(3, 1, 5, 2, 2),                            &
      subshell_term_LS(3, 2, 3, 2, 2),                            &
      subshell_term_LS(3, 3, 3, 2, 2),                            &
      subshell_term_LS(3, 1, 5,10, 2),                            &
      subshell_term_LS(3, 2, 3,10, 2),                            &
      subshell_term_LS(3, 3, 3,10, 2),                            &
      subshell_term_LS(3, 4, 3,10, 2),                            &
      subshell_term_LS(3, 4, 1, 2, 2),                            &
      subshell_term_LS(3, 5, 1,10, 2),                            &
      subshell_term_LS(3, 6, 1,10, 2),                            &
      subshell_term_LS(3, 5, 1, 2, 2),                            &
      subshell_term_LS(3, 6, 1, 2, 2),                            &
      subshell_term_LS(3, 7, 1,10, 2),                            &
      subshell_term_LS(3, 8, 1,10, 2),                            &
      subshell_term_LS(3, 9, 1,10, 2),                            &
      subshell_term_LS(3, 1, 3,12, 2),                            &
      subshell_term_LS(3, 2, 3,12, 2),                            &
      subshell_term_LS(3, 3, 1,12, 2),                            &
      subshell_term_LS(3, 3, 1,12, 2),                            &
      subshell_term_LS(3, 5, 1,12, 2),                            &
      subshell_term_LS(3, 6, 1,12, 2),                            &
      subshell_term_LS(3, 1, 3,14, 2),                            &
      subshell_term_LS(3, 2, 3,14, 2),                            &
      subshell_term_LS(3, 3, 1,14, 2),                            &
      subshell_term_LS(3, 4, 1,14, 2),                            &
      subshell_term_LS(3, 5, 1,14, 2),                            &
      subshell_term_LS(3, 6, 1,14, 2),                            &
      subshell_term_LS(3, 1, 3,16, 2),                            &
      subshell_term_LS(3, 2, 1,16, 2),                            &
      subshell_term_LS(3, 3, 1,16, 2),                            &
      subshell_term_LS(3, 1, 3,18, 2),                            &
      subshell_term_LS(3, 2, 1,18, 2),                            &
      subshell_term_LS(3, 3, 1,18, 2),                            &
      subshell_term_LS(3, 0, 1,20, 2),                            &
      subshell_term_LS(3, 0, 1,22, 2),                            &
      subshell_term_LS(3, 2, 1, 6, 0),                            &
      subshell_term_LS(3, 3, 1, 6, 0),                            &
      subshell_term_LS(3, 4, 1, 6, 0),                            &
      subshell_term_LS(3, 1, 5, 4, 0),                            &
      subshell_term_LS(3, 2, 3, 4, 0),                            &
      subshell_term_LS(3, 3, 3, 4, 0),                            &
      subshell_term_LS(3, 1, 3, 6, 0),                            &
      subshell_term_LS(3, 1, 5, 8, 0),                            &
      subshell_term_LS(3, 2, 3, 8, 0),                            &
      subshell_term_LS(3, 3, 3, 8, 0),                            &
      subshell_term_LS(3, 5, 1, 4, 0),                            &
      subshell_term_LS(3, 5, 1, 8, 0),                            &
      subshell_term_LS(3, 6, 1, 4, 0),                            &
      subshell_term_LS(3, 6, 1, 8, 0),                            &
      subshell_term_LS(3, 7, 1, 8, 0),                            &
      subshell_term_LS(3, 8, 1, 8, 0),                            &
      subshell_term_LS(3, 4, 3, 4, 0),                            &
      subshell_term_LS(3, 4, 3, 8, 0),                            &
      subshell_term_LS(3, 1, 3,10, 0),                            &
      subshell_term_LS(3, 2, 3,10, 0),                            &
      subshell_term_LS(3, 0, 1, 2, 0),                            &
      subshell_term_LS(3, 3, 1,10, 0),                            &
      subshell_term_LS(3, 4, 1,10, 0),                            &
      subshell_term_LS(3, 1, 7, 0, 0),                            &
      subshell_term_LS(3, 1, 5,12, 0),                            &
      subshell_term_LS(3, 2, 3, 0, 0),                            &
      subshell_term_LS(3, 2, 3,12, 0),                            &
      subshell_term_LS(3, 3, 3,12, 0),                            &
      subshell_term_LS(3, 3, 1, 0, 0),                            &
      subshell_term_LS(3, 4, 1,12, 0),                            &
      subshell_term_LS(3, 5, 1,12, 0),                            &
      subshell_term_LS(3, 4, 1, 0, 0),                            &
      subshell_term_LS(3, 6, 1,12, 0),                            &
      subshell_term_LS(3, 7, 1,12, 0),                            &
      subshell_term_LS(3, 1, 3,14, 0),                            &
      subshell_term_LS(3, 2, 1,14, 0),                            &
      subshell_term_LS(3, 3, 1,14, 0),                            &
      subshell_term_LS(3, 1, 3,16, 0),                            &
      subshell_term_LS(3, 2, 3,16, 0),                            &
      subshell_term_LS(3, 3, 1,16, 0),                            &
      subshell_term_LS(3, 4, 1,16, 0),                            &
      subshell_term_LS(3, 1, 1,18, 0),                            &
      subshell_term_LS(3, 2, 1,18, 0),                            &
      subshell_term_LS(3, 1, 3,20, 0),                            &
      subshell_term_LS(3, 2, 1,20, 0),                            &
      subshell_term_LS(3, 0, 1,24, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_g1 =(/                                      &
      subshell_term_LS(4, 0, 8, 8, 1)  /)
   type(subshell_term_LS), dimension(1:9), parameter ::   &
      term_LS_g2 =(/                                      &
      subshell_term_LS(4, 0, 9, 0, 0),                    &
      subshell_term_LS(4, 0, 7, 2, 2),                    &
      subshell_term_LS(4, 0, 7, 4, 0),                    &
      subshell_term_LS(4, 0, 7, 6, 2),                    &
      subshell_term_LS(4, 0, 7, 8, 0),                    &
      subshell_term_LS(4, 0, 7,10, 2),                    &
      subshell_term_LS(4, 0, 7,12, 0),                    &
      subshell_term_LS(4, 0, 7,14, 2),                    &
      subshell_term_LS(4, 0, 7,16, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_h1 =(/                                      &
      subshell_term_LS(5, 0,10,10, 1)  /)
   type(subshell_term_LS), dimension(1:11), parameter ::  &
      term_LS_h2 =(/                                      &
      subshell_term_LS(5, 0,11, 0, 0),                    &
      subshell_term_LS(5, 0, 9, 2, 2),                    &
      subshell_term_LS(5, 0, 9, 4, 0),                    &
      subshell_term_LS(5, 0, 9, 6, 2),                    &
      subshell_term_LS(5, 0, 9, 8, 0),                    &
      subshell_term_LS(5, 0, 9,10, 2),                    &
      subshell_term_LS(5, 0, 9,12, 0),                    &
      subshell_term_LS(5, 0, 9,14, 2),                    &
      subshell_term_LS(5, 0, 9,16, 0),                    &
      subshell_term_LS(5, 0, 9,18, 2),                    &
      subshell_term_LS(5, 0, 9,20, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_i1 =(/                                      &
      subshell_term_LS(6, 0,12,12, 1)  /)
   type(subshell_term_LS), dimension(1:13), parameter ::  &
      term_LS_i2 =(/                                      &
      subshell_term_LS(6, 0,13, 0, 0),                    &
      subshell_term_LS(6, 0,11, 2, 2),                    &
      subshell_term_LS(6, 0,11, 4, 0),                    &
      subshell_term_LS(6, 0,11, 6, 2),                    &
      subshell_term_LS(6, 0,11, 8, 0),                    &
      subshell_term_LS(6, 0,11,10, 2),                    &
      subshell_term_LS(6, 0,11,12, 0),                    &
      subshell_term_LS(6, 0,11,14, 2),                    &
      subshell_term_LS(6, 0,11,16, 0),                    &
      subshell_term_LS(6, 0,11,18, 2),                    &
      subshell_term_LS(6, 0,11,20, 0),                    &
      subshell_term_LS(6, 0,11,22, 2),                    &
      subshell_term_LS(6, 0,11,24, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_k1 =(/                                      &
      subshell_term_LS(7, 0,14,14, 1)  /)
   type(subshell_term_LS), dimension(1:15), parameter ::  &
      term_LS_k2 =(/                                      &
      subshell_term_LS(7, 0,15, 0, 0),                    &
      subshell_term_LS(7, 0,13, 2, 2),                    &
      subshell_term_LS(7, 0,13, 4, 0),                    &
      subshell_term_LS(7, 0,13, 6, 2),                    &
      subshell_term_LS(7, 0,13, 8, 0),                    &
      subshell_term_LS(7, 0,13,10, 2),                    &
      subshell_term_LS(7, 0,13,12, 0),                    &
      subshell_term_LS(7, 0,13,14, 2),                    &
      subshell_term_LS(7, 0,13,16, 0),                    &
      subshell_term_LS(7, 0,13,18, 2),                    &
      subshell_term_LS(7, 0,13,20, 0),                    &
      subshell_term_LS(7, 0,13,22, 2),                    &
      subshell_term_LS(7, 0,13,24, 0),                    &
      subshell_term_LS(7, 0,13,26, 2),                    &
      subshell_term_LS(7, 0,13,28, 0) /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_l1 =(/                                      &
      subshell_term_LS(8, 0,16,16, 1)  /)
   type(subshell_term_LS), dimension(1:17), parameter ::  &
      term_LS_l2 =(/                                      &
      subshell_term_LS(8, 0,17, 0, 0),                    &
      subshell_term_LS(8, 0,15, 2, 2),                    &
      subshell_term_LS(8, 0,15, 4, 0),                    &
      subshell_term_LS(8, 0,15, 6, 2),                    &
      subshell_term_LS(8, 0,15, 8, 0),                    &
      subshell_term_LS(8, 0,15,10, 2),                    &
      subshell_term_LS(8, 0,15,12, 0),                    &
      subshell_term_LS(8, 0,15,14, 2),                    &
      subshell_term_LS(8, 0,15,16, 0),                    &
      subshell_term_LS(8, 0,15,18, 2),                    &
      subshell_term_LS(8, 0,15,20, 0),                    &
      subshell_term_LS(8, 0,15,22, 2),                    &
      subshell_term_LS(8, 0,15,24, 0),                    &
      subshell_term_LS(8, 0,15,26, 2),                    &
      subshell_term_LS(8, 0,15,28, 0),                    &
      subshell_term_LS(8, 0,15,30, 2),                    &
      subshell_term_LS(8, 0,15,32, 0) /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_m1 =(/                                      &
      subshell_term_LS(9, 0,18,18, 1)  /)
   type(subshell_term_LS), dimension(1:19), parameter ::  &
      term_LS_m2 =(/                                      &
      subshell_term_LS(9, 0,19, 0, 0),                    &
      subshell_term_LS(9, 0,17, 2, 2),                    &
      subshell_term_LS(9, 0,17, 4, 0),                    &
      subshell_term_LS(9, 0,17, 6, 2),                    &
      subshell_term_LS(9, 0,17, 8, 0),                    &
      subshell_term_LS(9, 0,17,10, 2),                    &
      subshell_term_LS(9, 0,17,12, 0),                    &
      subshell_term_LS(9, 0,17,14, 2),                    &
      subshell_term_LS(9, 0,17,16, 0),                    &
      subshell_term_LS(9, 0,17,18, 2),                    &
      subshell_term_LS(9, 0,17,20, 0),                    &
      subshell_term_LS(9, 0,17,22, 2),                    &
      subshell_term_LS(9, 0,17,24, 0),                    &
      subshell_term_LS(9, 0,17,26, 2),                    &
      subshell_term_LS(9, 0,17,28, 0),                    &
      subshell_term_LS(9, 0,17,30, 2),                    &
      subshell_term_LS(9, 0,17,32, 0),                    &
      subshell_term_LS(9, 0,17,34, 2),                    &
      subshell_term_LS(9, 0,17,36, 0) /)
contains
   !
   !
   subroutine gettermjj(j_shell,N,jj,number)
   !--------------------------------------------------------------------
   !  This procedure return all allowed subshell terms
   !                                                 (2*j, Nr, 2*Q, 2*J)
   !  for given j^N which must be 1/2, 3/2, 5/2, 7/2 or 9/2.
   !
   !  The argument j_shell is 2*j
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)                                :: j_shell, N
      type(subshell_term), dimension(1:63), intent(out)  :: jj
      integer, intent(out)                               :: number
      !
      integer  :: M_Q, i, k
      !
      M_Q = N - (j_shell + 1)/2;  k = 0
      select case (j_shell)
      case (1)
         do i = 1,2
            if (mod(M_Q + terms_jj(i)%Q,2) == 0) then
               if (abs(M_Q) <= terms_jj(i)%Q) then
                  k = k + 1
                  jj(k)%j          = terms_jj(i)%j
                  jj(k)%Nr         = terms_jj(i)%Nr
                  jj(k)%Q          = terms_jj(i)%Q
                  jj(k)%subshellJ  = terms_jj(i)%subshellJ
                  jj(k)%nu         = (jj(k)%j+1)/2-jj(k)%Q
               end if
            end if
         end do
      case (3)
         do i = 3,5
            if (mod(M_Q + terms_jj(i)%Q,2) == 0) then
               if (abs(M_Q) <= terms_jj(i)%Q) then
                  k = k + 1
                  jj(k)%j          = terms_jj(i)%j
                  jj(k)%Nr         = terms_jj(i)%Nr
                  jj(k)%Q          = terms_jj(i)%Q
                  jj(k)%subshellJ  = terms_jj(i)%subshellJ
                  jj(k)%nu         = (jj(k)%j+1)/2-jj(k)%Q
               end if
            end if
         end do
      case (5)
         do i = 6,11
            if (mod(M_Q + terms_jj(i)%Q,2) == 0) then
               if (abs(M_Q) <= terms_jj(i)%Q) then
                  k = k + 1
                  jj(k)%j          = terms_jj(i)%j
                  jj(k)%Nr         = terms_jj(i)%Nr
                  jj(k)%Q          = terms_jj(i)%Q
                  jj(k)%subshellJ  = terms_jj(i)%subshellJ
                  jj(k)%nu         = (jj(k)%j+1)/2-jj(k)%Q
               end if
            end if
         end do
      case (7)
         do i = 12,25
            if (mod(M_Q + terms_jj(i)%Q,2) == 0) then
               if (abs(M_Q) <= terms_jj(i)%Q) then
                  k = k + 1
                  jj(k)%j          = terms_jj(i)%j
                  jj(k)%Nr         = terms_jj(i)%Nr
                  jj(k)%Q          = terms_jj(i)%Q
                  jj(k)%subshellJ  = terms_jj(i)%subshellJ
                  jj(k)%nu         = (jj(k)%j+1)/2-jj(k)%Q
               end if
            end if
         end do
      case (9)
         do i = 26,63
            if (mod(M_Q + terms_jj(i)%Q,2) == 0) then
               if (abs(M_Q) <= terms_jj(i)%Q) then
                  k = k + 1
                  jj(k)%j          = terms_jj(i)%j
                  jj(k)%Nr         = terms_jj(i)%Nr
                  jj(k)%Q          = terms_jj(i)%Q
                  jj(k)%subshellJ  = terms_jj(i)%subshellJ
                  jj(k)%nu         = (jj(k)%j+1)/2-jj(k)%Q
               end if
            end if
         end do
      case (11)
         select case (N)
         case (0)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_11_2(i)%j
            jj(k)%Nr        = terms_jj_11_2(i)%Nr
            jj(k)%Q         = terms_jj_11_2(i)%Q
            jj(k)%subshellJ = terms_jj_11_2(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (1)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_11_1(i)%j
            jj(k)%Nr        = terms_jj_11_1(i)%Nr
            jj(k)%Q         = terms_jj_11_1(i)%Q
            jj(k)%subshellJ = terms_jj_11_1(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (2)
            do i = 1,6
               k = k + 1
               jj(k)%j         = terms_jj_11_2(i)%j
               jj(k)%Nr        = terms_jj_11_2(i)%Nr
               jj(k)%Q         = terms_jj_11_2(i)%Q
               jj(k)%subshellJ = terms_jj_11_2(i)%subshellJ
               jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
            end do
         case default
            stop  "gettermLS(): program stop B."
         end select
      case (13)
         select case (N)
         case (0)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_13_2(i)%j
            jj(k)%Nr        = terms_jj_13_2(i)%Nr
            jj(k)%Q         = terms_jj_13_2(i)%Q
            jj(k)%subshellJ = terms_jj_13_2(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (1)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_13_1(i)%j
            jj(k)%Nr        = terms_jj_13_1(i)%Nr
            jj(k)%Q         = terms_jj_13_1(i)%Q
            jj(k)%subshellJ = terms_jj_13_1(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (2)
            do i = 1,7
               k = k + 1
               jj(k)%j         = terms_jj_13_2(i)%j
               jj(k)%Nr        = terms_jj_13_2(i)%Nr
               jj(k)%Q         = terms_jj_13_2(i)%Q
               jj(k)%subshellJ = terms_jj_13_2(i)%subshellJ
               jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
            end do
         case default
            stop  "gettermLS(): program stop C."
         end select
      case (15)
         select case (N)
         case (0)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_15_2(i)%j
            jj(k)%Nr        = terms_jj_15_2(i)%Nr
            jj(k)%Q         = terms_jj_15_2(i)%Q
            jj(k)%subshellJ = terms_jj_15_2(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (1)
            k = 1;       i = 1
            jj(k)%j         = terms_jj_15_1(i)%j
            jj(k)%Nr        = terms_jj_15_1(i)%Nr
            jj(k)%Q         = terms_jj_15_1(i)%Q
            jj(k)%subshellJ = terms_jj_15_1(i)%subshellJ
            jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
         case (2)
            do i = 1,8
               k = k + 1
               jj(k)%j         = terms_jj_15_2(i)%j
               jj(k)%Nr        = terms_jj_15_2(i)%Nr
               jj(k)%Q         = terms_jj_15_2(i)%Q
               jj(k)%subshellJ = terms_jj_15_2(i)%subshellJ
               jj(k)%nu        = (jj(k)%j+1)/2-jj(k)%Q
            end do
         case default
            stop  "gettermLS(): program stop D."
         end select
      case default
         print *, "j_shell = ", j_shell
         stop  "gettermjj(): program stop A."
      end select
      !
      number = k
      !
   end subroutine gettermjj
!
!***********************************************************************
!                                                                      *
      SUBROUTINE gettermLS(l_shell,N,LS,number)
!                                                                      *
!     This procedure returns all allowed subshell terms (l,w,Q,L,S)    *
!     for given l^N, N = 0, 1, 2 or 3.                                 *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!      USE jj2lsj_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)                                 :: l_shell, N
      type(subshell_term_LS), dimension(120), intent(out) :: LS
      integer, intent(out)                                :: number
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: M_Q, i, j
!-----------------------------------------------
      M_Q = N - 2* l_shell - 1;  j = 0
      select case (l_shell)
      case (0)
         do i = 1,2
            if (mod(M_Q + term_LS_s(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_s(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_s(i)%l_shell
                  LS(j)%w        = term_LS_s(i)%w
                  LS(j)%Q        = term_LS_s(i)%Q
                  LS(j)%LL       = term_LS_s(i)%LL
                  LS(j)%S        = term_LS_s(i)%S
               end if
            end if
         end do
      case (1)
         do i = 1,6
            if (mod(M_Q + term_LS_p(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_p(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_p(i)%l_shell
                  LS(j)%w        = term_LS_p(i)%w
                  LS(j)%Q        = term_LS_p(i)%Q
                  LS(j)%LL       = term_LS_p(i)%LL
                  LS(j)%S        = term_LS_p(i)%S
               end if
            end if
         end do
      case (2)
         do i = 1,32
            if (mod(M_Q + term_LS_d(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_d(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_d(i)%l_shell
                  LS(j)%w        = term_LS_d(i)%w
                  LS(j)%Q        = term_LS_d(i)%Q
                  LS(j)%LL       = term_LS_d(i)%LL
                  LS(j)%S        = term_LS_d(i)%S
               end if
            end if
         end do
      case (3)
         do i = 1,238
            if (mod(M_Q + term_LS_f(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_f(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_f(i)%l_shell
                  LS(j)%w        = term_LS_f(i)%w
                  LS(j)%Q        = term_LS_f(i)%Q
                  LS(j)%LL       = term_LS_f(i)%LL
                  LS(j)%S        = term_LS_f(i)%S
               end if
            end if
         end do
      case (4)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_g1(i)%l_shell
            LS(j)%w        = term_LS_g1(i)%w
            LS(j)%Q        = term_LS_g1(i)%Q
            LS(j)%LL       = term_LS_g1(i)%LL
            LS(j)%S        = term_LS_g1(i)%S
         case (2)
            do i = 1,9
               j = j + 1
               LS(j)%l_shell  = term_LS_g2(i)%l_shell
               LS(j)%w        = term_LS_g2(i)%w
               LS(j)%Q        = term_LS_g2(i)%Q
               LS(j)%LL       = term_LS_g2(i)%LL
               LS(j)%S        = term_LS_g2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop A."
         end select
      case (5)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_h1(i)%l_shell
            LS(j)%w        = term_LS_h1(i)%w
            LS(j)%Q        = term_LS_h1(i)%Q
            LS(j)%LL       = term_LS_h1(i)%LL
            LS(j)%S        = term_LS_h1(i)%S
         case (2)
            do i = 1,11
               j = j + 1
               LS(j)%l_shell  = term_LS_h2(i)%l_shell
               LS(j)%w        = term_LS_h2(i)%w
               LS(j)%Q        = term_LS_h2(i)%Q
               LS(j)%LL       = term_LS_h2(i)%LL
               LS(j)%S        = term_LS_h2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop B."
         end select
      case (6)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_i1(i)%l_shell
            LS(j)%w        = term_LS_i1(i)%w
            LS(j)%Q        = term_LS_i1(i)%Q
            LS(j)%LL       = term_LS_i1(i)%LL
            LS(j)%S        = term_LS_i1(i)%S
         case (2)
            do i = 1,13
               j = j + 1
               LS(j)%l_shell  = term_LS_i2(i)%l_shell
               LS(j)%w        = term_LS_i2(i)%w
               LS(j)%Q        = term_LS_i2(i)%Q
               LS(j)%LL       = term_LS_i2(i)%LL
               LS(j)%S        = term_LS_i2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop C."
         end select
      case (7)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_k1(i)%l_shell
            LS(j)%w        = term_LS_k1(i)%w
            LS(j)%Q        = term_LS_k1(i)%Q
            LS(j)%LL       = term_LS_k1(i)%LL
            LS(j)%S        = term_LS_k1(i)%S
         case (2)
            do i = 1,15
               j = j + 1
               LS(j)%l_shell  = term_LS_k2(i)%l_shell
               LS(j)%w        = term_LS_k2(i)%w
               LS(j)%Q        = term_LS_k2(i)%Q
               LS(j)%LL       = term_LS_k2(i)%LL
               LS(j)%S        = term_LS_k2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop D."
         end select
      case (8)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_l1(i)%l_shell
            LS(j)%w        = term_LS_l1(i)%w
            LS(j)%Q        = term_LS_l1(i)%Q
            LS(j)%LL       = term_LS_l1(i)%LL
            LS(j)%S        = term_LS_l1(i)%S
         case (2)
            do i = 1,17
               j = j + 1
               LS(j)%l_shell  = term_LS_l2(i)%l_shell
               LS(j)%w        = term_LS_l2(i)%w
               LS(j)%Q        = term_LS_l2(i)%Q
               LS(j)%LL       = term_LS_l2(i)%LL
               LS(j)%S        = term_LS_l2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop E."
         end select
      case (9)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_m1(i)%l_shell
            LS(j)%w        = term_LS_m1(i)%w
            LS(j)%Q        = term_LS_m1(i)%Q
            LS(j)%LL       = term_LS_m1(i)%LL
            LS(j)%S        = term_LS_m1(i)%S
         case (2)
            do i = 1,19
               j = j + 1
               LS(j)%l_shell  = term_LS_m2(i)%l_shell
               LS(j)%w        = term_LS_m2(i)%w
               LS(j)%Q        = term_LS_m2(i)%Q
               LS(j)%LL       = term_LS_m2(i)%LL
               LS(j)%S        = term_LS_m2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop F."
         end select
      case default
         stop  "gettermLS(): program stop G."
      end select
      number = j
      END SUBROUTINE gettermLS
!
!***********************************************************************
!                                                                      *
      FUNCTION gettermLSQ(l_shell,N,w,L,iS)          result (wa)
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                last update: February 2020   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_shell, N, w, L, iS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(subshell_term_LS), dimension(120) :: LS
      real (kind=dp) :: wa
      integer        :: number, i
!-----------------------------------------------
      CALL  gettermLS(l_shell,N,LS,number)
      do i = 1,number
         if (LS(i)%LL == L) then
            if (LS(i)%S == iS) then
               if (LS(i)%w == w) then
                  wa = LS(i)%Q
                  exit
               end if
            end if
         end if
      end do
      END FUNCTION gettermLSQ
!
      function coefLSjj(l_shell,N,w,Q,L,S,J,jm_shell,Nm,Qm,Jm,   &
                                               jp_shell,Qp,Jp) result(wa)
   !--------------------------------------------------------------------
   ! Returns the value of the LS-jj transformation matrix for a given
   ! set of quantum numbers.
   !
   ! Note that all (generalized) angular momentum quantum numbers except
   ! of l must be given twice the original numbers, i.e. for the quantum
   ! numbers Q, L, S, J, jm_shell, Qm, Jm, jp_shell, Qp, Jp.
   !
   ! Calls: coefLSjj2,  coefLSjjs.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: l_shell, N, w, Q, L, S, J,                &
                             jm_shell, Nm, Qm, Jm, jp_shell, Qp, Jp
      !
      real(kind=dp)       :: wa
      integer             :: NN, N1, N2, factor
      !
      wa = 0.0; factor = 1.0
      if (l_shell == 0)  then
         NN = N;   N1 = Nm;   N2 = N - Nm;
         if (S == J .and. jm_shell == 1 .and. NN == N1)  then
            wa = 1.0
         else if (S == J .and. jp_shell == 1 .and. NN == N2)  then
            wa = 1.0
         else
            wa = 0.0
         end if
      else if (jm_shell < jp_shell .or. jm_shell == jp_shell)  then
         if (l_shell > 0 .and. N > 2*l_shell +1)  then
            NN = 4*l_shell + 2 - N;   N1 = jm_shell + 1 - Nm
            N2 = jp_shell + 1 - N + Nm
            if (mod(2*l_shell+1-Q -((jm_shell+1)/2-Qm)         &
                 -((jp_shell+1)/2-Qp),4) /= 0) factor = - factor
         else
            NN = N;   N1 = Nm;   N2 = N - Nm;
         end if
         if (NN == 1  .or.  NN == 0)  then
            if (NN == 0 .and. N1 == 0 .and. N2 == 0) then
               wa = 1.0
            else if (N1 == 1 .and. N2 == 0 .and. J == Jm) then
               wa = 1.0
            else if (N1 == 0 .and. N2 == 1 .and. J == Jp) then
               wa = 1.0
            else
               wa = 0.0
            end if
         else if (NN == 2)  then
            if (J > Jm + Jp  .or.  J < abs(Jm - Jp)) then
               wa = 0.0
            else
               if (N1 == 2 .and. N2 == 0) then
                  wa = coefLSjj2 (l_shell,L,S,J,jm_shell,jm_shell)
               else if (N1 == 1 .and. N2 == 1) then
                  wa = coefLSjj2 (l_shell,L,S,J,jm_shell,jp_shell)
               else if (N1 == 0 .and. N2 == 2) then
                  wa = coefLSjj2 (l_shell,L,S,J,jp_shell,jp_shell)
               end if;
            end if
         else if (l_shell == 1  .or.  l_shell == 2  .or.  l_shell == 3)  then
            wa = coefLSjjs (l_shell,NN,w,Q,L,S,J,N1,Qm,Jm,Qp,Jp)
         end if
      else
         stop "coefLS(): program stop A."
      end if
      wa = factor * wa
      !
      end function coefLSjj
   !
   !
      function coefLSjj2 (l_shell,L,S,J,jm_shell,jp_shell)  result(wa)
   !--------------------------------------------------------------------
   ! Returns the value of the LS-jj transformation matrix
   !                                                (l^2 LSJ| j_1 j_2 J).
   !
   ! Note that all (generalized) angular momentum quantum numbers except
   ! of l must be given twice the original numbers, i.e. for the quantum
   ! numbers L, S, J, jm_shell, jp_shell.
   !
   ! Calls: wigner_9j_triangle, wigner_9j_symbol.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: l_shell, L, S, J, jm_shell, jp_shell
      real(kind=dp)       :: wa, RAC9
      !
      integer :: delta_J
      !
      wa = 0.0
      if (mod(L+S,4) /= 0)  then
         wa = 0.0
      elseif (mod(L+S+J,2) /= 0)  then
         wa = 0.0
      elseif (mod(jm_shell+jp_shell+J,2) /= 0)  then
         wa = 0.0
      else
         if (jm_shell == jp_shell) then
            if (mod(J,4) /= 0)  then
               wa = 0.0
            else
               wa = (jm_shell+1.0) * sqrt((L+1.0) * (S+1.0))
            end if
         else
            wa = sqrt(2 * (L+1.0) * (S+1.0) * (jm_shell+1.0) * (jp_shell+1.0))
         end if
         call nine(2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J,   &
                   1,delta_J,RAC9)
         if (delta_J /= 0)  then
            call nine(2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J,&
                      0,delta_J,RAC9)
            wa = wa * RAC9
!         delta_J = wigner_9j_triangle                                        &
!                      (2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J)
!         if (delta_J /= 0)  then
!            wa = wa * wigner_9j_symbol                                       &
!                      (2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J,.true.)

         end if
      end if
!
      end function coefLSjj2
!
      function coefLSjjs(lshell,N,w,Q,L,S,J,Nm,Qm,Jm,Qp,Jp) result(wa)
   !--------------------------------------------------------------------
   ! Returns the value of the LS-jj transformation matrix for a given
   ! set of quantum numbers.
   !
   ! Note that all (generalized) angular momentum quantum numbers except
   ! of l must be given twice the original numbers, i.e. for the quantum
   ! numbers Q, L, S, J, Qm, Jm, Qp, Jp.
   !
   ! Calls:
   !--------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: lshell, N, w, Q, L, S, J, Nm, Qm, Jm, Qp, Jp
      integer             :: i
      real(kind=dp)       :: wa
      !
      wa = 0.0
      if (lshell == 0) then
      else if (lshell == 1) then
         select case(N)
         case(3)
            ! Use data from the array LS_jj_p_3
            do  i = 1,LS_jj_number_p3
               if (w  == LS_jj_p_3(i)%w   .and. Q  == LS_jj_p_3(i)%Q  .and. &
                   L  == LS_jj_p_3(i)%L   .and. S  == LS_jj_p_3(i)%S  .and. &
                   J  == LS_jj_p_3(i)%J   .and. Nm == LS_jj_p_3(i)%Nm .and. &
                   Qm == LS_jj_p_3(i)%Qm  .and. Jm == LS_jj_p_3(i)%Jm .and. &
                   Qp == LS_jj_p_3(i)%Qp  .and. Jp == LS_jj_p_3(i)%Jp) then
                  wa = LS_jj_p_3(i)%factor * &
                       sqrt( 1.0*LS_jj_p_3(i)%nom/LS_jj_p_3(i)%denom )
                  return
               end if
            end do
         case(4)
            ! Use data from the array LS_jj_p_4
            do  i = 1,LS_jj_number_p4
               if (w  == LS_jj_p_4(i)%w   .and. Q  == LS_jj_p_4(i)%Q  .and. &
                   L  == LS_jj_p_4(i)%L   .and. S  == LS_jj_p_4(i)%S  .and. &
                   J  == LS_jj_p_4(i)%J   .and. Nm == LS_jj_p_4(i)%Nm .and. &
                   Qm == LS_jj_p_4(i)%Qm  .and. Jm == LS_jj_p_4(i)%Jm .and. &
                   Qp == LS_jj_p_4(i)%Qp  .and. Jp == LS_jj_p_4(i)%Jp) then
                  wa = LS_jj_p_4(i)%factor * &
                       sqrt( 1.0*LS_jj_p_4(i)%nom/LS_jj_p_4(i)%denom )
                  return
               end if
            end do
         case(5)
            ! Use data from the array LS_jj_p_5
            do  i = 1,LS_jj_number_p5
               if (w  == LS_jj_p_5(i)%w   .and. Q  == LS_jj_p_5(i)%Q  .and. &
                   L  == LS_jj_p_5(i)%L   .and. S  == LS_jj_p_5(i)%S  .and. &
                   J  == LS_jj_p_5(i)%J   .and. Nm == LS_jj_p_5(i)%Nm .and. &
                   Qm == LS_jj_p_5(i)%Qm  .and. Jm == LS_jj_p_5(i)%Jm .and. &
                   Qp == LS_jj_p_5(i)%Qp  .and. Jp == LS_jj_p_5(i)%Jp) then
                  wa = LS_jj_p_5(i)%factor * &
                       sqrt( 1.0*LS_jj_p_5(i)%nom/LS_jj_p_5(i)%denom )
                  return
               end if
            end do
         case(6)
            ! Use data from the array LS_jj_p_6
            do  i = 1,LS_jj_number_p6
               if (w  == LS_jj_p_6(i)%w   .and. Q  == LS_jj_p_6(i)%Q  .and. &
                   L  == LS_jj_p_6(i)%L   .and. S  == LS_jj_p_6(i)%S  .and. &
                   J  == LS_jj_p_6(i)%J   .and. Nm == LS_jj_p_6(i)%Nm .and. &
                   Qm == LS_jj_p_6(i)%Qm  .and. Jm == LS_jj_p_6(i)%Jm .and. &
                   Qp == LS_jj_p_6(i)%Qp  .and. Jp == LS_jj_p_6(i)%Jp) then
                  wa = LS_jj_p_6(i)%factor * &
                       sqrt( 1.0*LS_jj_p_6(i)%nom/LS_jj_p_6(i)%denom )
                  return
               end if
            end do
         case default
            stop "LS_jj_coefficient(): program stop A."
         end select
      else if (lshell == 2) then
         select case(N)
         case(3)
            ! Use data from the array LS_jj_d_3
            do  i = 1,LS_jj_number_d3
               if (w  == LS_jj_d_3(i)%w   .and. Q  == LS_jj_d_3(i)%Q  .and. &
                   L  == LS_jj_d_3(i)%L   .and. S  == LS_jj_d_3(i)%S  .and. &
                   J  == LS_jj_d_3(i)%J   .and. Nm == LS_jj_d_3(i)%Nm .and. &
                   Qm == LS_jj_d_3(i)%Qm  .and. Jm == LS_jj_d_3(i)%Jm .and. &
                   Qp == LS_jj_d_3(i)%Qp  .and. Jp == LS_jj_d_3(i)%Jp) then
                  wa = LS_jj_d_3(i)%factor * &
                       sqrt( 1.0*LS_jj_d_3(i)%nom/LS_jj_d_3(i)%denom )
                  return
               end if
            end do
         case(4)
            ! Use data from the array LS_jj_d_4
            do  i = 1,LS_jj_number_d4
               if (w  == LS_jj_d_4(i)%w   .and. Q  == LS_jj_d_4(i)%Q  .and. &
                   L  == LS_jj_d_4(i)%L   .and. S  == LS_jj_d_4(i)%S  .and. &
                   J  == LS_jj_d_4(i)%J   .and. Nm == LS_jj_d_4(i)%Nm .and. &
                   Qm == LS_jj_d_4(i)%Qm  .and. Jm == LS_jj_d_4(i)%Jm .and. &
                   Qp == LS_jj_d_4(i)%Qp  .and. Jp == LS_jj_d_4(i)%Jp) then
                  wa = LS_jj_d_4(i)%factor * &
                       sqrt( 1.0*LS_jj_d_4(i)%nom/LS_jj_d_4(i)%denom )
                  return
               end if
            end do
         case(5)
            ! Use data from the array LS_jj_d_5
            do  i = 1,LS_jj_number_d5
               if (w  == LS_jj_d_5(i)%w   .and. Q  == LS_jj_d_5(i)%Q  .and. &
                   L  == LS_jj_d_5(i)%L   .and. S  == LS_jj_d_5(i)%S  .and. &
                   J  == LS_jj_d_5(i)%J   .and. Nm == LS_jj_d_5(i)%Nm .and. &
                   Qm == LS_jj_d_5(i)%Qm  .and. Jm == LS_jj_d_5(i)%Jm .and. &
                   Qp == LS_jj_d_5(i)%Qp  .and. Jp == LS_jj_d_5(i)%Jp) then
                  wa = LS_jj_d_5(i)%factor * &
                       sqrt( 1.0*LS_jj_d_5(i)%nom/LS_jj_d_5(i)%denom )
                  return
               end if
            end do
         case(6)
            ! Use data from the array LS_jj_d_6
            do  i = 1,LS_jj_number_d6
               if (w  == LS_jj_d_6(i)%w   .and. Q  == LS_jj_d_6(i)%Q  .and. &
                   L  == LS_jj_d_6(i)%L   .and. S  == LS_jj_d_6(i)%S  .and. &
                   J  == LS_jj_d_6(i)%J   .and. Nm == LS_jj_d_6(i)%Nm .and. &
                   Qm == LS_jj_d_6(i)%Qm  .and. Jm == LS_jj_d_6(i)%Jm .and. &
                   Qp == LS_jj_d_6(i)%Qp  .and. Jp == LS_jj_d_6(i)%Jp) then
                  wa = LS_jj_d_6(i)%factor * &
                       sqrt( 1.0*LS_jj_d_6(i)%nom/LS_jj_d_6(i)%denom )
                  return
               end if
            end do
         case(7)
            ! Use data from the array LS_jj_d_7
            do  i = 1,LS_jj_number_d7
               if (w  == LS_jj_d_7(i)%w   .and. Q  == LS_jj_d_7(i)%Q  .and. &
                   L  == LS_jj_d_7(i)%L   .and. S  == LS_jj_d_7(i)%S  .and. &
                   J  == LS_jj_d_7(i)%J   .and. Nm == LS_jj_d_7(i)%Nm .and. &
                   Qm == LS_jj_d_7(i)%Qm  .and. Jm == LS_jj_d_7(i)%Jm .and. &
                   Qp == LS_jj_d_7(i)%Qp  .and. Jp == LS_jj_d_7(i)%Jp) then
                  wa = LS_jj_d_7(i)%factor * &
                       sqrt( 1.0*LS_jj_d_7(i)%nom/LS_jj_d_7(i)%denom )
                  return
               end if
            end do
         case(8)
            ! Use data from the array LS_jj_d_8
            do  i = 1,LS_jj_number_d8
               if (w  == LS_jj_d_8(i)%w   .and. Q  == LS_jj_d_8(i)%Q  .and. &
                   L  == LS_jj_d_8(i)%L   .and. S  == LS_jj_d_8(i)%S  .and. &
                   J  == LS_jj_d_8(i)%J   .and. Nm == LS_jj_d_8(i)%Nm .and. &
                   Qm == LS_jj_d_8(i)%Qm  .and. Jm == LS_jj_d_8(i)%Jm .and. &
                   Qp == LS_jj_d_8(i)%Qp  .and. Jp == LS_jj_d_8(i)%Jp) then
                  wa = LS_jj_d_8(i)%factor * &
                       sqrt( 1.0*LS_jj_d_8(i)%nom/LS_jj_d_8(i)%denom )
                  return
               end if
            end do
         case(9)
            ! Use data from the array LS_jj_d_9
            do  i = 1,LS_jj_number_d9
               if (w  == LS_jj_d_9(i)%w   .and. Q  == LS_jj_d_9(i)%Q  .and. &
                   L  == LS_jj_d_9(i)%L   .and. S  == LS_jj_d_9(i)%S  .and. &
                   J  == LS_jj_d_9(i)%J   .and. Nm == LS_jj_d_9(i)%Nm .and. &
                   Qm == LS_jj_d_9(i)%Qm  .and. Jm == LS_jj_d_9(i)%Jm .and. &
                   Qp == LS_jj_d_9(i)%Qp  .and. Jp == LS_jj_d_9(i)%Jp) then
                  wa = LS_jj_d_9(i)%factor * &
                       sqrt( 1.0*LS_jj_d_9(i)%nom/LS_jj_d_9(i)%denom )
                  return
               end if
            end do
         case(10)
            ! Use data from the array LS_jj_d_10
            do  i = 1,LS_jj_number_d10
               if (w  == LS_jj_d_10(i)%w   .and. Q  == LS_jj_d_10(i)%Q  .and. &
                   L  == LS_jj_d_10(i)%L   .and. S  == LS_jj_d_10(i)%S  .and. &
                   J  == LS_jj_d_10(i)%J   .and. Nm == LS_jj_d_10(i)%Nm .and. &
                   Qm == LS_jj_d_10(i)%Qm  .and. Jm == LS_jj_d_10(i)%Jm .and. &
                   Qp == LS_jj_d_10(i)%Qp  .and. Jp == LS_jj_d_10(i)%Jp) then
                  wa = LS_jj_d_10(i)%factor * &
                       sqrt( 1.0*LS_jj_d_10(i)%nom/LS_jj_d_10(i)%denom )
                  return
               end if
            end do
         case default
            stop "LS_jj_coefficient(): program stop B."
         end select
      else if (lshell == 3) then
         select case(N)
         case(3)
            ! Use data from the array LS_jj_f_3
            do  i = 1,LS_jj_number_f3
               if (w  == LS_jj_f_3(i)%w   .and. Q  == LS_jj_f_3(i)%Q  .and. &
                   L  == LS_jj_f_3(i)%L   .and. S  == LS_jj_f_3(i)%S  .and. &
                   J  == LS_jj_f_3(i)%J   .and. Nm == LS_jj_f_3(i)%Nm .and. &
                   Qm == LS_jj_f_3(i)%Qm  .and. Jm == LS_jj_f_3(i)%Jm .and. &
                   Qp == LS_jj_f_3(i)%Qp  .and. Jp == LS_jj_f_3(i)%Jp) then
                  wa = LS_jj_f_3(i)%factor * &
                       sqrt( 1.0*LS_jj_f_3(i)%nom/LS_jj_f_3(i)%denom )
                  return
               end if
            end do
         case(4)
            ! Use data from the array LS_jj_f_4
            do  i = 1,LS_jj_number_f4
               if (w  == LS_jj_f_4(i)%w   .and. Q  == LS_jj_f_4(i)%Q  .and. &
                   L  == LS_jj_f_4(i)%L   .and. S  == LS_jj_f_4(i)%S  .and. &
                   J  == LS_jj_f_4(i)%J   .and. Nm == LS_jj_f_4(i)%Nm .and. &
                   Qm == LS_jj_f_4(i)%Qm  .and. Jm == LS_jj_f_4(i)%Jm .and. &
                   Qp == LS_jj_f_4(i)%Qp  .and. Jp == LS_jj_f_4(i)%Jp) then
                  wa = LS_jj_f_4(i)%factor * &
                       sqrt( 1.0*LS_jj_f_4(i)%nom/LS_jj_f_4(i)%denom )
                  return
               end if
            end do
         case(5)
            ! Use data from the array LS_jj_f_5
            do  i = 1,LS_jj_number_f5
               if (w  == LS_jj_f_5(i)%w   .and. Q  == LS_jj_f_5(i)%Q  .and. &
                   L  == LS_jj_f_5(i)%L   .and. S  == LS_jj_f_5(i)%S  .and. &
                   J  == LS_jj_f_5(i)%J   .and. Nm == LS_jj_f_5(i)%Nm .and. &
                   Qm == LS_jj_f_5(i)%Qm  .and. Jm == LS_jj_f_5(i)%Jm .and. &
                   Qp == LS_jj_f_5(i)%Qp  .and. Jp == LS_jj_f_5(i)%Jp) then
                  wa = LS_jj_f_5(i)%factor * &
                       sqrt( 1.0*LS_jj_f_5(i)%nom/LS_jj_f_5(i)%denom )
                  return
               end if
            end do
         case(6)
            ! Use data from the array LS_jj_f_6
            do  i = 1,LS_jj_number_f6
               if (w  == LS_jj_f_6(i)%w   .and. Q  == LS_jj_f_6(i)%Q  .and. &
                   L  == LS_jj_f_6(i)%L   .and. S  == LS_jj_f_6(i)%S  .and. &
                   J  == LS_jj_f_6(i)%J   .and. Nm == LS_jj_f_6(i)%Nm .and. &
                   Qm == LS_jj_f_6(i)%Qm  .and. Jm == LS_jj_f_6(i)%Jm .and. &
                   Qp == LS_jj_f_6(i)%Qp  .and. Jp == LS_jj_f_6(i)%Jp) then
                  wa = LS_jj_f_6(i)%factor * &
                       sqrt( 1.0*LS_jj_f_6(i)%nom/LS_jj_f_6(i)%denom )
                  return
               end if
            end do
         case(7)
            ! Use data from the array LS_jj_f_7
            do  i = 1, LS_jj_number_f7, 1
               if (w == LS_jj_f_7(i)%w) then
                  if(Q == LS_jj_f_7(i)%Q) then
                     if(L == LS_jj_f_7(i)%L) then
                        if(S == LS_jj_f_7(i)%S) then
                           if(J == LS_jj_f_7(i)%J) then
                              if(Nm == LS_jj_f_7(i)%Nm) then
                                 if(Qm == LS_jj_f_7(i)%Qm .and. &
                                    Jm == LS_jj_f_7(i)%Jm .and. &
                                    Qp == LS_jj_f_7(i)%Qp .and. &
                                    Jp == LS_jj_f_7(i)%Jp) then
                                       wa = 1.0 * LS_jj_f_7(i)%factor  * &
                                            sqrt(1.0*LS_jj_f_7(i)%nom)/ &
                                            sqrt(1.0*LS_jj_f_7(i)%denom )
                                       return
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
            end do
         case default
            stop "coefLSjjs(): program stop C."
         end select
      else
         stop "coefLSjjs(): program stop D."
      end if
      !
      !
      end function coefLSjjs
!
!     -------------------------------------------------------------
!      J T H N
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION JTHN(K,N,I)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, N, K
!-----------------------------------------------
      IF(N.EQ.1) THEN
        JTHN=MOD(K,I)
      ELSEIF(N.EQ.2) THEN
        JTHN=MOD(K/I,I)
      ELSEIF(N.EQ.3) THEN
        JTHN=MOD(K/(I*I),I)
      ELSEIF(N.EQ.4) THEN
        JTHN=MOD(K/(I*I*I),I)
      ELSEIF(N.EQ.5) THEN
        JTHN=MOD(K/(I*I*I*I),I)
      ELSE
        WRITE(6,'(A)') ' ERROR IN JTHN '
        STOP
      ENDIF
      RETURN
      END FUNCTION JTHN
!
!     -------------------------------------------------------------
!      g e t J M i n M a x
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!
      SUBROUTINE getJMinMax(l,N,J_Min,J_Max)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, N
      integer, intent(out):: J_Min,J_Max
!-----------------------------------------------
!
      if(N==0 .or. 4*l+2-N==0) then
         J_Min = 0
         J_Max = 0
      else if(N==1 .or. 4*l+2-N==1) then
         J_Min = iabs(2*l-1)
         J_Max = 2*l+1
      else if(N==2 .or. 4*l+2-N==2) then
         J_Min = 0
         J_Max = 4*l
      else if(l == 1) then
         if(N==3) then
            J_Min = 1
            J_Max = 5
         else
            WRITE(6,'(A)') 'ERROR IN getJMinMax'
            WRITE(6,'(A6,2I5)') 'l,N = ',l,N
            STOP
         end if
      else if(l==2) then
!GG!         if(N==3) then
         if(N==3 .or. N==7) then
            J_Min = 1
            J_Max = 11
!GG!         else if(N==4) then
         else if(N==4 .or. N==6) then
            J_Min = 0
            J_Max = 12
         else if(N==5) then
            J_Min = 1
            J_Max = 13
         else
            WRITE(6,'(A)') 'ERROR IN getJMinMax'
            WRITE(6,'(A6,2I5)') 'l,N = ',l,N
            STOP
         end if
      else if(l==3) then
!GG!         if(N==3) then
         if(N==3 .or. N==11) then
            J_Min = 1
            J_Max = 17
!GG!         else if(N==4) then
         else if(N==4 .or. N==10) then
            J_Min = 0
            J_Max = 20
!GG!         else if(N==5) then
         else if(N==5 .or. N==9) then
            J_Min = 1
            J_Max = 23
!GG!         else if(N==6) then
         else if(N==6 .or. N==8) then
            J_Min = 0
            J_Max = 24
         else if(N==7) then
            J_Min = 1
            J_Max = 25
         else
            WRITE(6,'(A)') 'ERROR IN getJMinMax'
            WRITE(6,'(A6,2I5)') 'l,N = ',l,N
            STOP
         end if
      else
         WRITE(6,'(A)') 'ERROR IN getJMinMax'
         WRITE(6,'(A6,2I5)') 'l,N = ',l,N
         STOP
      end if
      RETURN
      END SUBROUTINE getJMinMax
!
!     -------------------------------------------------------------
!      g e t j j J M i n M a x
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!
      SUBROUTINE getjjJMinMax(j,N,J_Min,J_Max)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: j, N
      integer, intent(out):: J_Min,J_Max
!-----------------------------------------------
!
      if(N==0 .or. j+1-N==0) then
         J_Min = 0
         J_Max = 0
      else if(N==1 .or. j+1-N==1) then
         J_Min = j
         J_Max = j
      else if(N==2 .or. j+1-N==2) then
         J_Min = 0
         J_Max = 2*(j-1)
      else if(j == 5) then
         if(N==3) then
!GG!!            J_Min = 5
            J_Min = 3
            J_Max = 9
         else
            WRITE(6,'(A)') 'ERROR IN getjjJMinMax'
            WRITE(6,'(A8,2I5)') '2*j,N = ',j,N
            STOP
         end if
      else if(j==7) then
         if(N==3 .or. N==5) then
!GG!!            J_Min = 7
            J_Min = 3
            J_Max = 15
         else if(N==4) then
            J_Min = 0
            J_Max = 16
         else
            WRITE(6,'(A)') 'ERROR IN getjjJMinMax'
            WRITE(6,'(A8,2I5)') '2*j,N = ',j,N
            STOP
         end if
      else if(j==9) then
!GG!         if(N==3) then
         if(N==3 .or. N==7) then
            J_Min = 3
            J_Max = 21
!GG!         else if(N==4) then
         else if(N==4 .or. N==6) then
            J_Min = 0
            J_Max = 24
         else if(N==5) then
            J_Min = 3
            J_Max = 25
         else
            WRITE(6,'(A)') 'ERROR IN getjjJMinMax'
            WRITE(6,'(A8,2I5)') '2*j,N = ',j,N
            STOP
         end if
      else
         WRITE(6,'(A)') 'ERROR IN getjjJMinMax'
         WRITE(6,'(A8,2I5)') '2*j,N = ',j,N
         STOP
      end if
      RETURN
      END SUBROUTINE getjjJMinMax

end module Coupling_inside_shell
