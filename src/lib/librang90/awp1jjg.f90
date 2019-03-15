!*******************************************************************
!                                                                  *
      SUBROUTINE AWP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,AW)
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
      USE CONS_C,          ONLY: ZERO, TWO, TENTH, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjj2_I
      USE ixjtik_I
      USE izas1_I
      USE rumtjj_I
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
      REAL(DOUBLE), INTENT(OUT)              :: AW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK1, KK2, IE, IQ, IQ2, IQ3, IQM, IT, ITP, ITG, IBTT
      INTEGER,      DIMENSION(7) :: IBT
      REAL(DOUBLE)               :: ENQP, D1, SI1, W
      REAL(DOUBLE), DIMENSION(3) :: BT
!-----------------------------------------------
      AW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=BK2+BK2+TENTH*BK2
      IQ2=QM2*TWO+QM2*TENTH
      IQ3=QM3*TWO+QM3*TENTH
      IQ=IQ2+IQ2
      IF(ID(3) > 37) RETURN
      IF(ITJJ2(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ) == 0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM <= IBT(7)) THEN
          IF(IXJTIK(IK(3),KK1,KK2,ID(6),IK(6),IBT(6)) /= 0) THEN
            IBT(1)=IT
            BT(2)=DBLE(IBT(6))/TWO
            BT(1)=DBLE(IBT(7))/TWO
            CALL A1JJ(IK,BK,IBT,BT,QM1,D1)
            IF(DABS(D1) > EPS) THEN
              CALL W1JJG(K1,QM2,QM3,IBT,BT,ID,BD,W)
              IF(DABS(W) > EPS) THEN
                D1=D1*W
                CALL SIXJ(IK(3),KK1,KK2,ID(6),IK(6),IBT(6),0,SI1)
                ENQP=ENQP+D1*SI1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      END DO
      AW=ENQP*DSQRT(DBLE((KK2+1)))
      IE=KK2+IK(6)+ID(6)
      IF(((IE/4)*4) /= IE)AW=-AW
      RETURN
      END SUBROUTINE AWP1JJG
