!*******************************************************************
!                                                                  *
      SUBROUTINE ONESCALAR2INT(JJA,JJB,JA,JB,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,ITREXG,IXJTIK,PERKO2,              *
!                        RECO,RECO2                                *
!                                                                  *
!   Written by  G. Gaigalas                   NIST, December 2015  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE m_C
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoonescalar_I
      USE perko2_I
      USE reco2_I
!      USE gg12_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IAT,JAA,JBB,NN,IB1,II,INN,IFAZ
      REAL(DOUBLE) :: QM1,QM2,REC,WW
!-----------------------------------------------
      INTERACT=0
      IF(JA.EQ.JB) RETURN
      IF(JA.LT.JB) THEN
        JAA=JA
        JBB=JB
      ELSE
        JAA=JB
        JBB=JA
      END IF
      CALL RECOONESCALAR(-1,JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT.EQ.0)RETURN
      QM1=HALF
      QM2=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      IF(ID1(3).NE.ID2(3)) RETURN
      CALL RECO2(JAA,JBB,ID2(3),0,IAT,REC)
      IF(IAT.EQ.0)RETURN
!GG INT
      INTERACT=1
!
!      CALL GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)
!      IF(DABS(WW).GT.EPS) THEN
!         CALL RECO2(JAA,JBB,ID2(3),1,IAT,REC)
!         COEFF=WW*REC*DSQRT(DBLE(ID1(3)+1))
!         NN=0
!         IB1=JBB-1
!         DO II=JAA,IB1
!           INN=JLIST(II)
!           NN=NQ1(INN)+NN
!         END DO
!         IF((NN/2)*2.EQ.NN) COEFF=-COEFF
!iGGGGG         IF(JA.GT.JB) COEFF=-COEFF
!         COEFF=-COEFF
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
!        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
!        IF((IFAZ/4)*4.NE.IFAZ)COEFF=-COEFF
!      END IF
      RETURN
      END SUBROUTINE ONESCALAR2INT
