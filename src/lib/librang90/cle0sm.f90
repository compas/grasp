!*******************************************************************
!                                                                  *
      SUBROUTINE CLE0SM(Q,QM,S,C,CM,A)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   S  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                            December 1993   *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, TENTH, ONE, TWO, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      USE c1e0sm_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, S, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IIS
!-----------------------------------------------
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IIS=TWO*S+TENTH
      IF(ITTK(IIQ,IIC,IIS).EQ.0)RETURN
      IF(S.LT.EPS) THEN
       IF((Q+TENTH).LT.DABS(QM))RETURN
        IF((C+TENTH).LT.DABS(CM))RETURN
        IF(DABS(Q-C).GT.EPS)RETURN
        IF(DABS(QM-CM).GT.EPS)RETURN
        A=ONE
      ELSE
        CALL C1E0SM(Q,QM,C,CM,A)
      END IF
      RETURN
      END SUBROUTINE CLE0SM
