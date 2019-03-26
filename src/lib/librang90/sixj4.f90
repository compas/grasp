!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ4(JC,JE,JD,JB,JF,ITIK,SI)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | JC/2  JE/2  JD/2 |                                         *
!     | JB/2  JF/2    4  |                                         *
!                                                                  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, HALF, ONE, TWO, THREE, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      USE dracah_I
      USE sixj2_I
      USE sixj3_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER             :: JC, JE, JD, JB, JF
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE)        :: SI
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: A, C, E, D, B, F, X1, X2, X3, S2, S3
!-----------------------------------------------
      SI = ZERO
      IF (ITIK /= 0) THEN
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(JC,JE,JD,JB,JF,8) == 0) RETURN
      ENDIF
      IF (IXJTIK(JC,JE,JD,JB,JF,6) == 0) THEN
         CALL DRACAH (JC, JE, JF, JB, JD, 8, SI)
         IF (MOD(JC + JE + JF + JB,4) /= 0) SI = -SI
      ELSE
         A = THREE
         C = DBLE(JC)*HALF
         E = DBLE(JE)*HALF
         D = DBLE(JD)*HALF
         B = DBLE(JB)*HALF
         F = DBLE(JF)*HALF
         X1 = A*DSQRT((A+B+E+TWO)*(A-B+E+ONE)*(A+B-E+ONE)*((-         &
            A)+B+E)*(A+C+F+TWO)*(A-C+F+ONE)*(A+C-F+ONE)*(             &
            (-A)+C+F))
         X2 = (A+ONE)*DSQRT((A+B+E+ONE)*(A-B+E)*(A+B-E)*((-A)         &
             + B+E+ONE)*(A+C+F+ONE)*(A-C+F)*(A+C-F)*((-A)+C           &
             + F+ONE))
         X3 = (TWO*A+ONE)*(TWO*(A*(A+ONE)*D*(D+ONE)-B*(B+ONE)*C*(C+   &
            ONE)-E*(E+ONE)*F*(F+ONE))+(A*(A+ONE)-B*(B+ONE)-E*(E       &
             +ONE))*(A*(A+ONE)-C*(C+ONE)-F*(F+ONE)))
         IF (DABS(X2) < EPS) THEN
            S2 = ZERO
         ELSE
            CALL SIXJ2 (JC, JE, JD, JB, JF, 0, S2)
         ENDIF
         IF (DABS(X3) < EPS) THEN
            S3 = ZERO
         ELSE
            CALL SIXJ3 (JC, JE, JD, JB, JF, 0, S3)
         ENDIF
         SI = (X3*S3 - X2*S2)/X1
      ENDIF
      RETURN
      END SUBROUTINE SIXJ4
