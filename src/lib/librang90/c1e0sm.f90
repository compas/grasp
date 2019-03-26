!*******************************************************************
!                                                                  *
      SUBROUTINE C1E0SM(Q,QM,C,CM,A)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, TENTH, ONE, TWO, THREE, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IS, IG
!-----------------------------------------------
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,2) == 0)RETURN
      IF(DABS(QM-CM) > EPS) RETURN
      IF((Q+TENTH) < DABS(QM)) RETURN
      IF((C+TENTH) < DABS(CM)) RETURN
      IF(DABS(QM) <= EPS) THEN
       IS=Q+C+ONE+TENTH
       IF((IS/2)*2 /= IS) RETURN
      END IF
      IG=Q-C+TWO+TENTH
      IF(IG <= 0) RETURN
      IF(IG > 3) RETURN
      IF (IG == 1) THEN
        A=DSQRT(((C+CM)*(C-CM))/((TWO*C-ONE)*C))
      ELSE IF (IG == 2) THEN
        A=CM/DSQRT(C*(C+ONE))
      ELSE IF (IG == 3) THEN
        A=-DSQRT(((C+CM+ONE)*(C-CM+ONE))/((C+ONE)*(TWO*C+THREE)))
      END IF
      RETURN
      END SUBROUTINE C1E0SM
