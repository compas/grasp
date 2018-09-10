!
!***********************************************************************
!                                                                      *
      MODULE jj2lsj_C 
!                                                                      *
!     This module contains the (numerical) values of the subshell      *
!     terms in LS-coupling and in it are define some global variables  *
!     for jj2LSJ program.                                              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  LONG, DOUBLE 
!
! Define some global data for the jj2LSJ transformation program
!
      type :: subshell_term_LS
         integer :: l_shell ! angular momentum l
         integer :: w       ! w
         integer :: Q       ! Subshell total quasispin 2*Q
         integer :: LL      ! Subshell total angular momentum 2*L
         integer :: S       ! Subshell Total sinp 2*S
      end type subshell_term_LS
!
      type :: LS_jj_me
!         sequence
         integer       :: w, Q, L, S, J, Nm, Qm, Jm, Qp, Jp
         integer       :: factor
         integer(LONG) :: nom, denom
      end type LS_jj_me
!
! Define the values of all LS-jj transformation coefficients
!
      integer, parameter :: LS_jj_number_p3 =  10, LS_jj_number_p4  =    9, &
                            LS_jj_number_p5 =   2, LS_jj_number_p6  =    1, &
                            LS_jj_number_d3 =  65, LS_jj_number_d4  =  166, &
                            LS_jj_number_d5 = 184, LS_jj_number_d6  =  166, &
                            LS_jj_number_d7 =  65, LS_jj_number_d8  =   19, &
                            LS_jj_number_d9 =   2, LS_jj_number_d10 =    1, &
                            LS_jj_number_f3 = 216, LS_jj_number_f4  = 1210, &
                            LS_jj_number_f5 =3799, LS_jj_number_f6  = 7313, &
                            LS_jj_number_f7 =8003
!
! Define all possible antisymetric subshell terms in LS-coupling
!
   type(subshell_term_LS), dimension(1:2), parameter ::   &
      term_LS_s =(/                                       &
      subshell_term_LS(0, 0, 0, 0, 1),                    &
      subshell_term_LS(0, 0, 1, 0, 0)  /)
   type(subshell_term_LS), dimension(1:6), parameter ::   &
      term_LS_p =(/                                       &
      subshell_term_LS(1, 0, 0, 0, 3),                    &
      subshell_term_LS(1, 0, 2, 2, 1),                    &
      subshell_term_LS(1, 0, 0, 4, 1),                    &
      subshell_term_LS(1, 0, 3, 0, 0),                    &
      subshell_term_LS(1, 0, 1, 2, 2),                    &
      subshell_term_LS(1, 0, 1, 4, 0)  /)
   type(subshell_term_LS), dimension(1:32), parameter ::  &
      term_LS_d =(/                                       &
      subshell_term_LS(2, 0, 0, 0, 5),                    &
      subshell_term_LS(2, 0, 0, 0, 1),                    &
      subshell_term_LS(2, 0, 2, 2, 3),                    &
      subshell_term_LS(2, 0, 2, 2, 1),                    &
      subshell_term_LS(2, 0, 4, 4, 1),                    &
      subshell_term_LS(2, 0, 2, 4, 1),                    &
      subshell_term_LS(2, 0, 0, 4, 3),                    &
      subshell_term_LS(2, 0, 0, 4, 1),                    &
      subshell_term_LS(2, 0, 2, 6, 3),                    &
      subshell_term_LS(2, 0, 2, 6, 1),                    &
      subshell_term_LS(2, 0, 0, 6, 1),                    &
      subshell_term_LS(2, 0, 2, 8, 1),                    &
      subshell_term_LS(2, 0, 0, 8, 3),                    &
      subshell_term_LS(2, 0, 0, 8, 1),                    &
      subshell_term_LS(2, 0, 2,10, 1),                    &
      subshell_term_LS(2, 0, 0,12, 1),                    &
      subshell_term_LS(2, 0, 5, 0, 0),                    &
      subshell_term_LS(2, 0, 1, 0, 0),                    &
      subshell_term_LS(2, 0, 3, 2, 2),                    &
      subshell_term_LS(2, 0, 1, 2, 2),                    &
      subshell_term_LS(2, 0, 1, 4, 4),                    &
      subshell_term_LS(2, 0, 1, 4, 2),                    &
      subshell_term_LS(2, 0, 3, 4, 0),                    &
      subshell_term_LS(2, 0, 1, 4, 0),                    &
      subshell_term_LS(2, 0, 3, 6, 2),                    &
      subshell_term_LS(2, 0, 1, 6, 2),                    &
      subshell_term_LS(2, 0, 1, 6, 0),                    &
      subshell_term_LS(2, 0, 1, 8, 2),                    &
      subshell_term_LS(2, 0, 3, 8, 0),                    &
      subshell_term_LS(2, 0, 1, 8, 0),                    &
      subshell_term_LS(2, 0, 1,10, 2),                    &
      subshell_term_LS(2, 0, 1,12, 0)  /)
   type(subshell_term_LS), dimension(1:238), parameter :: &
      term_LS_f =(/                                       &
      subshell_term_LS(3, 0, 0, 0, 7),                    &
      subshell_term_LS(3, 0, 2, 2, 5),                    &
      subshell_term_LS(3, 0, 0, 4, 5),                    &
      subshell_term_LS(3, 0, 2, 6, 5),                    &
      subshell_term_LS(3, 0, 0, 8, 5),                    &
      subshell_term_LS(3, 0, 2,10, 5),                    &
      subshell_term_LS(3, 0, 0,12, 5),                    &
      subshell_term_LS(3, 1, 4, 0, 3),                    &
      subshell_term_LS(3, 2, 0, 0, 3),                    &
      subshell_term_LS(3, 1, 2, 2, 3),                    &
      subshell_term_LS(3, 2, 2, 2, 3),                    &
      subshell_term_LS(3, 1, 4, 4, 3),                    &
      subshell_term_LS(3, 2, 2, 4, 3),                    &
      subshell_term_LS(3, 3, 2, 4, 3),                    &
      subshell_term_LS(3, 4, 0, 4, 3),                    &
      subshell_term_LS(3, 5, 0, 4, 3),                    &
      subshell_term_LS(3, 6, 0, 4, 3),                    &
      subshell_term_LS(3, 1, 4, 6, 3),                    &
      subshell_term_LS(3, 2, 2, 6, 3),                    &
      subshell_term_LS(3, 3, 2, 6, 3),                    &
      subshell_term_LS(3, 4, 2, 6, 3),                    &
      subshell_term_LS(3, 5, 0, 6, 3),                    &
      subshell_term_LS(3, 1, 4, 8, 3),                    &
      subshell_term_LS(3, 2, 2, 8, 3),                    &
      subshell_term_LS(3, 3, 2, 8, 3),                    &
      subshell_term_LS(3, 4, 2, 8, 3),                    &
      subshell_term_LS(3, 5, 0, 8, 3),                    &
      subshell_term_LS(3, 6, 0, 8, 3),                    &
      subshell_term_LS(3, 7, 0, 8, 3),                    &
      subshell_term_LS(3, 1, 2,10, 3),                    &
      subshell_term_LS(3, 2, 2,10, 3),                    &
      subshell_term_LS(3, 3, 2,10, 3),                    &
      subshell_term_LS(3, 4, 0,10, 3),                    &
      subshell_term_LS(3, 5, 0,10, 3),                    &
      subshell_term_LS(3, 1, 4,12, 3),                    &
      subshell_term_LS(3, 2, 2,12, 3),                    &
      subshell_term_LS(3, 3, 2,12, 3),                    &
      subshell_term_LS(3, 4, 0,12, 3),                    &
      subshell_term_LS(3, 5, 0,12, 3),                    &
      subshell_term_LS(3, 1, 2,14, 3),                    &
      subshell_term_LS(3, 2, 2,14, 3),                    &
      subshell_term_LS(3, 3, 0,14, 3),                    &
      subshell_term_LS(3, 1, 2,16, 3),                    &
      subshell_term_LS(3, 2, 0,16, 3),                    &
      subshell_term_LS(3, 3, 0,16, 3),                    &
      subshell_term_LS(3, 0, 2,18, 3),                    &
      subshell_term_LS(3, 0, 0,20, 3),                    &
      subshell_term_LS(3, 1, 0, 0, 1),                    &
      subshell_term_LS(3, 2, 0, 0, 1),                    &
      subshell_term_LS(3, 1, 4, 2, 1),                    &
      subshell_term_LS(3, 2, 2, 2, 1),                    &
      subshell_term_LS(3, 3, 2, 2, 1),                    &
      subshell_term_LS(3, 4, 2, 2, 1),                    &
      subshell_term_LS(3, 5, 0, 2, 1),                    &
      subshell_term_LS(3, 1, 4, 4, 1),                    &
      subshell_term_LS(3, 2, 4, 4, 1),                    &
      subshell_term_LS(3, 3, 2, 4, 1),                    &
      subshell_term_LS(3, 4, 2, 4, 1),                    &
      subshell_term_LS(3, 5, 2, 4, 1),                    &
      subshell_term_LS(3, 6, 0, 4, 1),                    &
      subshell_term_LS(3, 7, 0, 4, 1),                    &
      subshell_term_LS(3, 1, 6, 6, 1),                    &
      subshell_term_LS(3, 2, 4, 6, 1),                    &
      subshell_term_LS(3, 3, 2, 6, 1),                    &
      subshell_term_LS(3, 4, 2, 6, 1),                    &
      subshell_term_LS(3, 5, 2, 6, 1),                    &
      subshell_term_LS(3, 6, 2, 6, 1),                    &
      subshell_term_LS(3, 7, 2, 6, 1),                    &
      subshell_term_LS(3, 8, 0, 6, 1),                    &
      subshell_term_LS(3, 9, 0, 6, 1),                    &
      subshell_term_LS(3,10, 0, 6, 1),                    &
      subshell_term_LS(3, 1, 4, 8, 1),                    &
      subshell_term_LS(3, 2, 4, 8, 1),                    &
      subshell_term_LS(3, 3, 2, 8, 1),                    &
      subshell_term_LS(3, 4, 2, 8, 1),                    &
      subshell_term_LS(3, 5, 2, 8, 1),                    &
      subshell_term_LS(3, 6, 2, 8, 1),                    &
      subshell_term_LS(3, 7, 0, 8, 1),                    &
      subshell_term_LS(3, 8, 0, 8, 1),                    &
      subshell_term_LS(3, 9, 0, 8, 1),                    &
      subshell_term_LS(3,10, 0, 8, 1),                    &
      subshell_term_LS(3, 1, 4,10, 1),                    &
      subshell_term_LS(3, 2, 4,10, 1),                    &
      subshell_term_LS(3, 3, 2,10, 1),                    &
      subshell_term_LS(3, 4, 2,10, 1),                    &
      subshell_term_LS(3, 5, 2,10, 1),                    &
      subshell_term_LS(3, 6, 2,10, 1),                    &
      subshell_term_LS(3, 7, 2,10, 1),                    &
      subshell_term_LS(3, 8, 0,10, 1),                    &
      subshell_term_LS(3, 9, 0,10, 1),                    &
      subshell_term_LS(3, 1, 4,12, 1),                    &
      subshell_term_LS(3, 2, 2,12, 1),                    &
      subshell_term_LS(3, 3, 2,12, 1),                    &
      subshell_term_LS(3, 4, 2,12, 1),                    &
      subshell_term_LS(3, 5, 2,12, 1),                    &
      subshell_term_LS(3, 6, 0,12, 1),                    &
      subshell_term_LS(3, 7, 0,12, 1),                    &
      subshell_term_LS(3, 8, 0,12, 1),                    &
      subshell_term_LS(3, 9, 0,12, 1),                    &
      subshell_term_LS(3, 1, 4,14, 1),                    &
      subshell_term_LS(3, 2, 2,14, 1),                    &
      subshell_term_LS(3, 3, 2,14, 1),                    &
      subshell_term_LS(3, 4, 2,14, 1),                    &
      subshell_term_LS(3, 5, 2,14, 1),                    &
      subshell_term_LS(3, 6, 0,14, 1),                    &
      subshell_term_LS(3, 7, 0,14, 1),                    &
      subshell_term_LS(3, 1, 4,16, 1),                    &
      subshell_term_LS(3, 2, 2,16, 1),                    &
      subshell_term_LS(3, 3, 2,16, 1),                    &
      subshell_term_LS(3, 4, 0,16, 1),                    &
      subshell_term_LS(3, 5, 0,16, 1),                    &
      subshell_term_LS(3, 1, 2,18, 1),                    &
      subshell_term_LS(3, 2, 2,18, 1),                    &
      subshell_term_LS(3, 3, 0,18, 1),                    &
      subshell_term_LS(3, 4, 0,18, 1),                    &
      subshell_term_LS(3, 1, 2,20, 1),                    &
      subshell_term_LS(3, 2, 0,20, 1),                    &
      subshell_term_LS(3, 0, 2,22, 1),                    &
      subshell_term_LS(3, 0, 0,24, 1),                    &
      subshell_term_LS(3, 0, 1, 6, 6),                    &
      subshell_term_LS(3, 1, 3, 4, 4),                    &
      subshell_term_LS(3, 2, 1, 4, 4),                    &
      subshell_term_LS(3, 3, 1, 4, 4),                    &
      subshell_term_LS(3, 1, 3, 6, 4),                    &
      subshell_term_LS(3, 2, 1, 6, 4),                    &
      subshell_term_LS(3, 1, 3, 8, 4),                    &
      subshell_term_LS(3, 2, 1, 8, 4),                    &
      subshell_term_LS(3, 3, 1, 8, 4),                    &
      subshell_term_LS(3, 0, 1, 2, 4),                    &
      subshell_term_LS(3, 1, 1,10, 4),                    &
      subshell_term_LS(3, 2, 1,10, 4),                    &
      subshell_term_LS(3, 0, 3, 0, 4),                    &
      subshell_term_LS(3, 1, 3,12, 4),                    &
      subshell_term_LS(3, 2, 1,12, 4),                    &
      subshell_term_LS(3, 0, 1,14, 4),                    &
      subshell_term_LS(3, 0, 1,16, 4),                    &
      subshell_term_LS(3, 1, 5, 6, 2),                    &
      subshell_term_LS(3, 2, 3, 6, 2),                    &
      subshell_term_LS(3, 6, 1, 6, 2),                    &
      subshell_term_LS(3, 8, 1, 6, 2),                    &
      subshell_term_LS(3, 1, 3, 4, 2),                    &
      subshell_term_LS(3, 2, 3, 4, 2),                    &
      subshell_term_LS(3, 3, 1, 4, 2),                    &
      subshell_term_LS(3, 4, 1, 4, 2),                    &
      subshell_term_LS(3, 3, 3, 6, 2),                    &
      subshell_term_LS(3, 5, 1, 6, 2),                    &
      subshell_term_LS(3, 1, 3, 8, 2),                    &
      subshell_term_LS(3, 2, 3, 8, 2),                    &
      subshell_term_LS(3, 4, 1, 8, 2),                    &
      subshell_term_LS(3, 5, 1, 8, 2),                    &
      subshell_term_LS(3, 5, 1, 4, 2),                    &
      subshell_term_LS(3, 4, 3, 6, 2),                    &
      subshell_term_LS(3, 7, 1, 6, 2),                    &
      subshell_term_LS(3, 9, 1, 6, 2),                    &
      subshell_term_LS(3, 3, 3, 8, 2),                    &
      subshell_term_LS(3, 6, 1, 8, 2),                    &
      subshell_term_LS(3, 7, 1, 8, 2),                    &
      subshell_term_LS(3, 1, 5, 2, 2),                    &
      subshell_term_LS(3, 2, 3, 2, 2),                    &
      subshell_term_LS(3, 3, 3, 2, 2),                    &
      subshell_term_LS(3, 1, 5,10, 2),                    &
      subshell_term_LS(3, 2, 3,10, 2),                    &
      subshell_term_LS(3, 3, 3,10, 2),                    &
      subshell_term_LS(3, 4, 3,10, 2),                    &
      subshell_term_LS(3, 4, 1, 2, 2),                    &
      subshell_term_LS(3, 5, 1,10, 2),                    &
      subshell_term_LS(3, 6, 1,10, 2),                    &
      subshell_term_LS(3, 5, 1, 2, 2),                    &
      subshell_term_LS(3, 6, 1, 2, 2),                    &
      subshell_term_LS(3, 7, 1,10, 2),                    &
      subshell_term_LS(3, 8, 1,10, 2),                    &
      subshell_term_LS(3, 9, 1,10, 2),                    &
      subshell_term_LS(3, 1, 3,12, 2),                    &
      subshell_term_LS(3, 2, 3,12, 2),                    &
      subshell_term_LS(3, 3, 1,12, 2),                    &
      subshell_term_LS(3, 3, 1,12, 2),                    &
      subshell_term_LS(3, 5, 1,12, 2),                    &
      subshell_term_LS(3, 6, 1,12, 2),                    &
      subshell_term_LS(3, 1, 3,14, 2),                    &
      subshell_term_LS(3, 2, 3,14, 2),                    &
      subshell_term_LS(3, 3, 1,14, 2),                    &
      subshell_term_LS(3, 4, 1,14, 2),                    &
      subshell_term_LS(3, 5, 1,14, 2),                    &
      subshell_term_LS(3, 6, 1,14, 2),                    &
      subshell_term_LS(3, 1, 3,16, 2),                    &
      subshell_term_LS(3, 2, 1,16, 2),                    &
      subshell_term_LS(3, 3, 1,16, 2),                    &
      subshell_term_LS(3, 1, 3,18, 2),                    &
      subshell_term_LS(3, 2, 1,18, 2),                    &
      subshell_term_LS(3, 3, 1,18, 2),                    &
      subshell_term_LS(3, 0, 1,20, 2),                    &
      subshell_term_LS(3, 0, 1,22, 2),                    &
      subshell_term_LS(3, 2, 1, 6, 0),                    &
      subshell_term_LS(3, 3, 1, 6, 0),                    &
      subshell_term_LS(3, 4, 1, 6, 0),                    &
      subshell_term_LS(3, 1, 5, 4, 0),                    &
      subshell_term_LS(3, 2, 3, 4, 0),                    &
      subshell_term_LS(3, 3, 3, 4, 0),                    &
      subshell_term_LS(3, 1, 3, 6, 0),                    &
      subshell_term_LS(3, 1, 5, 8, 0),                    &
      subshell_term_LS(3, 2, 3, 8, 0),                    &
      subshell_term_LS(3, 3, 3, 8, 0),                    &
      subshell_term_LS(3, 5, 1, 4, 0),                    &
      subshell_term_LS(3, 5, 1, 8, 0),                    &
      subshell_term_LS(3, 6, 1, 4, 0),                    &
      subshell_term_LS(3, 6, 1, 8, 0),                    &
      subshell_term_LS(3, 7, 1, 8, 0),                    &
      subshell_term_LS(3, 8, 1, 8, 0),                    &
      subshell_term_LS(3, 4, 3, 4, 0),                    &
      subshell_term_LS(3, 4, 3, 8, 0),                    &
      subshell_term_LS(3, 1, 3,10, 0),                    &
      subshell_term_LS(3, 2, 3,10, 0),                    &
      subshell_term_LS(3, 0, 1, 2, 0),                    &
      subshell_term_LS(3, 3, 1,10, 0),                    &
      subshell_term_LS(3, 4, 1,10, 0),                    &
      subshell_term_LS(3, 1, 7, 0, 0),                    &
      subshell_term_LS(3, 1, 5,12, 0),                    &
      subshell_term_LS(3, 2, 3, 0, 0),                    &
      subshell_term_LS(3, 2, 3,12, 0),                    &
      subshell_term_LS(3, 3, 3,12, 0),                    &
      subshell_term_LS(3, 3, 1, 0, 0),                    &
      subshell_term_LS(3, 4, 1,12, 0),                    &
      subshell_term_LS(3, 5, 1,12, 0),                    &
      subshell_term_LS(3, 4, 1, 0, 0),                    &
      subshell_term_LS(3, 6, 1,12, 0),                    &
      subshell_term_LS(3, 7, 1,12, 0),                    &
      subshell_term_LS(3, 1, 3,14, 0),                    &
      subshell_term_LS(3, 2, 1,14, 0),                    &
      subshell_term_LS(3, 3, 1,14, 0),                    &
      subshell_term_LS(3, 1, 3,16, 0),                    &
      subshell_term_LS(3, 2, 3,16, 0),                    &
      subshell_term_LS(3, 3, 1,16, 0),                    &
      subshell_term_LS(3, 4, 1,16, 0),                    &
      subshell_term_LS(3, 1, 1,18, 0),                    &
      subshell_term_LS(3, 2, 1,18, 0),                    &
      subshell_term_LS(3, 1, 3,20, 0),                    &
      subshell_term_LS(3, 2, 1,20, 0),                    &
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
      END MODULE jj2lsj_C 
