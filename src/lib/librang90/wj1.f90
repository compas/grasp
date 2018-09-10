!*******************************************************************
!                                                                  *
      SUBROUTINE WJ1(IK,BK,ID,BD,K2,QM1,QM2,WJ)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 24  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                      N       (k2)  N'     +-                     *
!     ELEMENT:       (j  QJ:: W   ::j  QJ)  -+                     *
!                                           ++                     *
!                                           -- S5(1.47),(1.48),    *
!                                                (1.49),(1.50).    *
!                                                                  *
!     SUBROUTINE CALLED: CLE0SM,C1E1SM,RWJJ                        *
!                                                                  *
!   Written by  G. Gaigalas                                        *
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
      USE mes_I
      USE w1jjg_I
      USE cle0sm_I
      USE c1e1sm_I
      USE rwjj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K2
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WJ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: K1, IQ, IQM2
      REAL(DOUBLE) :: A, QQ, W, WK1
!-----------------------------------------------
      WJ=ZERO
      IF(ID(3) == 9) THEN
        IF(MAX0(IK(4),ID(4)) < 3) THEN
          IF(IK(1) < 300) CALL MES(56)
          IF(ID(1) < 300) CALL MES(56)
          IQM2=QM2+QM2+QM2*EPS
          IF((ID(4)+IQM2) > 2) CALL MES(2)
          CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
          RETURN
        ELSE
          PRINT*, "ERROR in  WJ1"
          STOP
        ENDIF
      ELSEIF(ID(3) > 9) THEN
        IQM2=QM2+QM2+QM2*EPS
        IF((ID(4)+IQM2) > 2) CALL MES(2)
        CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
        RETURN
      ENDIF
      QQ=QM1+QM2
      IF(DABS(QQ) >= EPS) THEN
        IF(((K2/2)*2) /= K2)RETURN
        IQ=QQ+QQ*TENTH
        IF((IK(4)-ID(4)-2*IQ) /= 0)RETURN
        CALL C1E1SM(BD(1),BD(3),QQ,BK(1),BK(3),A)
        IF(DABS(A) < EPS)RETURN
        CALL RWJJ(IK(3),IK(1),ID(1),1,K2,W)
        WJ=A*W/DSQRT(TWO*BK(1)+ONE)
      ELSE IF(IK(4) /= ID(4))THEN
        RETURN
      ELSE IF(K2 == 0) THEN
        IF(ID(1) /= IK(1))RETURN
        IF(QM1 < EPS) THEN
          A=DBLE(ID(3)+1-ID(4))
        ELSE
          A=-DBLE(ID(4))
        END IF
        WJ=A*DSQRT(DBLE(IK(6)+1)/DBLE(IK(3)+1))
      ELSE
        K1=1
        IF(((K2/2)*2) /= K2)K1=0
        WK1=DBLE(K1)
        CALL CLE0SM(BD(1),BD(3),WK1,BK(1),BK(3),A)
        IF(DABS(A) < EPS)RETURN
        CALL RWJJ(IK(3),IK(1),ID(1),K1,K2,W)
        A=A*W
        WJ=A/DSQRT(TWO*TWO*BK(1)+TWO)
        IF(QM1 >= EPS)RETURN
        IF(((K2/2)*2) /= K2)WJ=-WJ
      END IF
      RETURN
      END SUBROUTINE WJ1
