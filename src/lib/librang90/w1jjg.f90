!*******************************************************************
!                                                                  *
      SUBROUTINE W1JJG(K1,QM1,QM2,IK,BK,ID,BD,WW)
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
      USE CONS_C,          ONLY: ZERO, TENTH, HALF, TWO, EPS
      USE ribojj_C
      USE ribojj9_C
      USE ribojj11_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE rumtjj_I
      USE ixjtik_I
      USE ittk_I
      USE a1jj_I
      USE sixj_I
      USE izas1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK1,IQ,IQM,IE,IT,ITK,ITD,ITP,ITP1,ITG,ITG1,IBTT
      INTEGER,      DIMENSION(7) :: IBT
      REAL(DOUBLE)               :: ENQP, D1, SI1, W
      REAL(DOUBLE), DIMENSION(3) :: BT
!-----------------------------------------------
      WW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK1=K1*2
      IQ=QM2*TWO+QM2*TENTH
      IF(ID(3) > 37) RETURN
      IF(ITTK(ID(6),IK(6),KK1) == 0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3) == 9) THEN
        IF(ID(4) > 2) CALL MES(1)
        IF(IK(4) > 2) CALL MES(1)
        ITK=ITK-300
        ITD=ITD-300
        ITP1=IMPNJJ9(ITK)
        ITP=IMPNJJ9(ITD)
        IF(ITP1 /= ITP)RETURN
        ITG1=IMGNJJ9(ITK)
        ITG=IMGNJJ9(ITD)
       ELSE
        IF(ID(4) > 2) CALL MES(1)
        IF(IK(4) > 2) CALL MES(1)
        ITP1=IMPNJJ11(ITK)
        ITP=IMPNJJ11(ITD)
        IF(ITP1 /= ITP)RETURN
        ITG1=IMGNJJ11(ITK)
        ITG=IMGNJJ11(ITD)
      ENDIF
      IF(ITG1 /= ITG)RETURN
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM <= IBT(7)) THEN
          IF(IXJTIK(IK(3),IK(3),KK1,ID(6),IK(6),IBT(6)) /= 0) THEN
            IBT(1)=IT
            BT(2)=DBLE(IBT(6))/TWO
            BT(1)=DBLE(IBT(7))/TWO
            CALL A1JJ(IK,BK,IBT,BT,QM1,D1)
            IF(DABS(D1) > EPS) THEN
              CALL A1JJ(IBT,BT,ID,BD,QM2,W)
              IF(DABS(W) > EPS) THEN
                D1=D1*W
                CALL SIXJ(IK(3),IK(3),KK1,ID(6),IK(6),IBT(6),0,SI1)
                ENQP=ENQP+D1*SI1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      END DO
      WW=ENQP*DSQRT(DBLE((KK1+1)))
      IE=KK1+IK(6)+ID(6)
      IF(((IE/4)*4) /= IE)WW=-WW
      RETURN
      END SUBROUTINE W1JJG
