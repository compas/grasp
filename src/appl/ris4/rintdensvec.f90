!***********************************************************************
!                                                                      *
      SUBROUTINE RINTDENSVEC (I,J,DINT1VEC,NRNUC)
!                                                                      *
!   The value of DINT1VEC is an approximation to:                      *
!                                                                      *
!                                                                      *
!      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) |                     *
!                       I     J       I     J                          *
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by Per Jonsson  and Jörgen Ekman                           *
!   Last revision: 24 Dec 2013                                         *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE def_C,            ONLY: pi
      USE grid_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,N), INTENT(OUT) :: DINT1VEC
      INTEGER  :: I, J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L, NRNUC
!-----------------------------------------------
!
      DO L = 2, NRNUC
         DINT1VEC(I,J,L) = (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))/          &
                   (4.0D00*PI*R(L)*R(L))
      END DO
!
      RETURN
      END SUBROUTINE RINTDENSVEC
