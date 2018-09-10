!*******************************************************************
!                                                                  *
      SUBROUTINE C1E1SM(Q,QM,SM,C,CM,A)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  1  CM  I  *
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
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, SM, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                    :: IE, IIQ, IIC
      REAL(DOUBLE), DIMENSION(2) :: GC
!-----------------------------------------------
      GC(1)=ONE
      GC(2)=-ONE
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
      IF(DABS(QM+SM-CM).GT.EPS)RETURN
      IF((Q+TENTH).LT.DABS(QM))RETURN
      IF((C+TENTH).LT.DABS(CM))RETURN
      IE=0
      IF(DABS(SM-ONE).LT.EPS)IE=1
      IF(DABS(SM+ONE).LT.EPS)IE=2
      IF(IE.EQ.0)RETURN
      IF(DABS(Q+ONE-C).LT.EPS) THEN
        A=DSQRT((C+GC(IE)*CM-ONE)*(C+GC(IE)*CM)/ &
        ((TWO*C-ONE)*TWO*C))
      ELSE IF(DABS(Q-C).LT.EPS) THEN
        A=-GC(IE)*DSQRT((C+GC(IE)*CM)*(C-GC(IE)*CM+ONE)/ &
        ((C+ONE)*TWO*C))
      ELSE IF(DABS(Q-ONE-C).GT.EPS) THEN
        RETURN
      ELSE
        A=DSQRT((C-GC(IE)*CM+ONE)*(C-GC(IE)*CM+TWO)/ &
        ((TWO*C+TWO)*(TWO*C+THREE)))
      END IF
      RETURN
      END SUBROUTINE C1E1SM
