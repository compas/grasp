!*******************************************************************
!                                                                  *
      SUBROUTINE C0T5S(Q,QM,SM,C,CM,A)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                ---          ---  *
!                                                I  Q  1/2  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:              I              I  *
!                                                I  QM  SM  CM  I  *
!                                                ---          ---  *
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
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C,          ONLY: ZERO, TENTH , HALF, ONE, TWO, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, SM, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
!      DIMENSION GC(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                    :: IIQ, IIC, IE
      REAL(DOUBLE), DIMENSION(2) :: GC
!-----------------------------------------------
      GC(1)=ONE
      GC(2)=-ONE
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,1) == 0)RETURN
      IF(DABS(QM+SM-CM) > EPS)RETURN
      IF((HALF+TENTH) < DABS(SM))RETURN
      IF((Q+TENTH) < DABS(QM))RETURN
      IF((C+TENTH) < DABS(CM))RETURN
      IE=DABS(HALF-SM)+ONE+TENTH
      IF(DABS(Q+HALF-C) < EPS) THEN
        A=DSQRT((C+GC(IE)*CM)/(TWO*C))
      ELSE
        IF(DABS(Q-HALF-C) > EPS)RETURN
        A=-GC(IE)*DSQRT((C-GC(IE)*CM+ONE)/(TWO*C+TWO))
      ENDIF
      RETURN
      END SUBROUTINE C0T5S
