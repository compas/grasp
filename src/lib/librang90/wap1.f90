!*******************************************************************
!                                                                  *
      SUBROUTINE WAP1(IK,BK,ID,BD,K1,BK2,QM1,QM2,QM3,WA)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 23  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                    N       (k1)   (j) (k2)   N'    +-            *
!     ELEMENT:     (j QJ::[ W    * A   ]    ::j QJ)  -+            *
!                                                    ++            *
!                                                    -- B17 (2.3)  *
!                                                                  *
!     SUBROUTINE CALLED: C0T5S,ITJJ3,IXJTIK,IZAS1,RUMTJJ,SIXJ,     *
!                        RMEAJJ,WJ1                                *
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
      USE mes_I
      USE wap1jjg_I
      USE izas1_I
      USE itjj3_I
      USE rumtjj_I
      USE ixjtik_I
      USE c0T5S_I
      USE rmeajj_I
      USE wj1_I
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
      REAL(DOUBLE)               :: ENQP, D1, S, SI, W
      REAL(DOUBLE), DIMENSION(3) :: BT
!-----------------------------------------------
      WA=ZERO
      IF(ID(3) == 9) THEN
        IF(MAX0(IK(4),ID(4)) < 3) THEN
          IF(IK(1) < 300) CALL MES(55)
          IF(ID(1) < 300) CALL MES(55)
          CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
          RETURN
        ENDIF
      ELSEIF(ID(3) > 9) THEN
        CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
        RETURN
      ENDIF
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK1=K1*2
      IQ3=QM3*TWO
      KK2=BK2+BK2+TENTH*BK2
      IF(ITJJ3(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM > IBT(7)) CYCLE
        IF(IXJTIK(KK1,IK(3),KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL C0T5S(BD(1),BD(3),QM3,BT(1),BT(3),D1)
        IF(DABS(D1) < EPS) CYCLE
        CALL RMEAJJ(ID(3),IBT(1),IBT(7),IBT(6),ID(1),ID(7),ID(6),S)
        IF(DABS(S) < EPS) CYCLE
        D1=D1*S
        CALL WJ1(IK,BK,IBT,BT,K1,QM1,QM2,W)
        IF(DABS(W) < EPS) CYCLE
        D1=D1*W
        CALL SIXJ(KK1,IK(3),KK2,ID(6),IK(6),IBT(6),0,SI)
        D1=D1*SI/DSQRT(DBLE(IBT(7)+1))
        ENQP=ENQP+D1
      END DO
      WA=ENQP*DSQRT(DBLE(KK2+1))
      IE=KK2+IK(6)+ID(6)+2
      IF(((IE/4)*4) /= IE)WA=-WA
      RETURN
      END SUBROUTINE WAP1
