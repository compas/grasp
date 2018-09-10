!*******************************************************************
!                                                                  *
      MODULE ribojj9_C 
!                                                                  *
!   This module is need for librang90.                             *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
! 
      IMPLICIT NONE
      INTEGER, DIMENSION(6) :: IMPTJJ9, IMGTJJ9, IMPNJJ9, IMGNJJ9
      DATA IMPTJJ9/301,5*302/
      DATA IMGTJJ9/301,5*306/
      DATA IMPNJJ9/302,5*301/
      DATA IMGNJJ9/306,5*301/
      END MODULE ribojj9_C 
