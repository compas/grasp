!*******************************************************************
!                                                                  *
      SUBROUTINE WW1(IK,BK,ID,BD,K2,QM1,QM2,QM3,QM4,WW)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 25  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                    N      (k2)   (k2) (0)   N'     +-            *
!     ELEMENT      (j QJ::[W   *  W    ]   ::j QJ)   -+            *
!                                                    ++            *
!                                                    -- B17 (2.4)  *
!                                                                  *
!     SUBROUTINE CALLED: ITJJ,IXJTIK,IZAS1,RUMTJJ,WJ1              *
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
      USE itjj_I
      USE rumtjj_I
      USE ixjtik_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K2
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK2,KK6,IQ,IQM,IQ3,IQ4,IE,IE1,IT,ITP,ITG,IBTT
      INTEGER, DIMENSION(7)      :: IBT
      REAL(DOUBLE)               :: ENQP, D1, W
      REAL(DOUBLE), DIMENSION(3) :: BT
!-----------------------------------------------
      WW=ZERO
      IF(ID(6) /= IK(6))RETURN
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK2=K2*2
      IQ3=QM3*TWO+QM3*TENTH
      IQ4=QM4*TWO+QM4*TENTH
      IQ=IQ3+IQ4
      IF(ITJJ(IK,ID,0,BK,BD,IBT,BT,KK6,ITP,ITG,IQ) == 0)RETURN
      IE1=KK2-IK(6)
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM > IBT(7)) CYCLE
        IF(IXJTIK(KK2,KK2,0,ID(6),IK(6),IBT(6)) == 0) CYCLE
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL WJ1(IK,BK,IBT,BT,K2,QM1,QM2,D1)
        IF(DABS(D1) < EPS) CYCLE
        CALL WJ1(IBT,BT,ID,BD,K2,QM3,QM4,W)
        IF(DABS(W) < EPS) CYCLE
        D1=D1*W
        IE=IE1+IBT(6)
        IF(((IE/4)*4) /= IE)D1=-D1
        ENQP=ENQP+D1
      END DO
      WW=ENQP/DSQRT(DBLE(KK2+1)*(IK(6)+1))
      RETURN
      END SUBROUTINE WW1
