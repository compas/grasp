!*******************************************************************
!                                                                  *
      SUBROUTINE WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
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
      USE CONS_C,          ONLY: ZERO, TENTH, TWO, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE izas1_I
      USE itjj3_I
      USE rumtjj_I
      USE ixjtik_I
      USE a1jj_I
      USE w1jjg_I
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: BK2, QM1, QM2, QM3
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK1, KK2, IE, IQ3, IQM, IT, ITP, ITG, IBTT
      INTEGER,      DIMENSION(7) :: IBT
      REAL(DOUBLE)               :: ENQP, D1, SI1, W
      REAL(DOUBLE), DIMENSION(3) :: BT
!-----------------------------------------------
      WA=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=BK2+BK2+TENTH*BK2
      IQ3=QM3*TWO+QM3*TENTH
      IF(ID(3) > 37) RETURN
      IF(ITJJ3(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM > IBT(7)) CYCLE
        IF(IXJTIK(KK1,IK(3),KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL A1JJ(IBT,BT,ID,BD,QM3,D1)
        IF(DABS(D1) < EPS) CYCLE
        CALL W1JJG(K1,QM1,QM2,IK,BK,IBT,BT,W)
        IF(DABS(W) < EPS) CYCLE
        D1=D1*W
        CALL SIXJ(KK1,IK(3),KK2,ID(6),IK(6),IBT(6),0,SI1)
        ENQP=ENQP+D1*SI1
        END DO
      WA=ENQP*DSQRT(DBLE((KK2+1)))
      IE=KK2+IK(6)+ID(6)
      IF(((IE/4)*4) /= IE)WA=-WA
      RETURN
      END SUBROUTINE WAP1JJG
