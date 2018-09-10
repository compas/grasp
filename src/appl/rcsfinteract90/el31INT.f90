!*******************************************************************
!                                                                  *
      SUBROUTINE EL31INT(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI, &
                                                           INTERACT)
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
      USE orb_C
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reco_I
      USE reco2_I
      USE perko2_I
      USE itrexg_I
      USE ixjtik_I
      USE coulom_I
      USE snrc_I
!      USE sixj_I
!      USE speak_I
!      USE gg1222_I
!      USE talk_I
!      USE cxk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,II,I2,I3,IAT,IIA,IIB,IIC,IID,IKK,IP1,IG1, &
!                  IBRD,IBRE,IFAZ,IFAZFRCS,INN,JAA,JBB,JB1,J12,L1, &
                  IBRD,IBRE,IFAZ,INN,JAA,JBB,JB1,J12,L1, &
                  L2,KRA,ND1,ND2,NE1,NE2,N,NN
      INTEGER :: INTERACT1, INTERACT2
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
!      REAL(DOUBLE)          :: AA,AB,A1,BB,SI,QM1,QM2,QM3,QM4,RECC
!      REAL(DOUBLE)          :: AA,AB,A1,BB,RECC
      REAL(DOUBLE)          :: A1,RECC
!      REAL(DOUBLE), DIMENSION(12) :: S
!-----------------------------------------------
      INTERACT = 0
      IF(NPEEL.LE.1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(JA.GT.JB) THEN
        JAA=JB
        JBB=JA
      ELSE
        JAA=JA
        JBB=JB
      END IF
      CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
!      QM1=HALF
!      QM2=HALF
!      QM3=-HALF
!      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      CALL RECO2(JAA,JBB,J(2),0,IAT,RECC)
      IF(IAT.EQ.0)RETURN
      IP1=ITREXG(J(1),J(1),J(1),J(2),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
!      CALL RECO2(JAA,JBB,J(2),1,IAT,RECC)
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=IABS(KAPS(1))
        KS(2)=IABS(KAPS(2))
        KS(3)=IABS(KAPS(3))
        KS(4)=IABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD .LE. 0)RETURN
      END IF
! * * *                      * * *                      * * *
!     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +
!           1211                                1112
!
      DO I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI .EQ. 1) THEN
          INTERACT1 = 0
          CALL COULOM(L2,L1,L1,L1,ID2(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(DABS(A1).LT.EPS) CYCLE
          INTERACT1 = 1
!          A1=-A1
        END IF
!        AB=ZERO
        INTERACT2 = 0
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF((IFAZ/2)*2.NE.IFAZ) CYCLE
          IF(IXJTIK(J(2),J(1),KRA*2,J(1),J(1),J12*2).EQ.0)CYCLE
          INTERACT2 = INTERACT2 + 1
!          CALL GG1222(IK2,IK1,BK2,BK1,ID2,ID1,BD2,BD1,J12,  &
!                      QM1,QM2,QM3,QM4,AA)
!          IF(DABS(AA).LT.EPS) CYCLE
!          CALL SIXJ(J(2),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
!          AA=AA*SI*DSQRT(DBLE(I3))
!          IFAZ=2*J(1)+KRA*2+J12*2
!          IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
!          AB=AB+AA
        END DO
        IF(INTERACT2 .EQ. 0) CYCLE
!        AB=AB*RECC
!        IF(DABS(AB).LT.EPS) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
!        IFAZFRCS = 1
!        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
!        IF((IFAZ/4)*4.NE.IFAZ)IFAZFRCS=-IFAZFRCS
!
!        NN=0
!        JB1=JBB-1
!        DO II=JAA,JB1
!          INN=JLIST(II)
!          NN=NQ1(INN)+NN
!        END DO
!        IF((NN/2)*2.EQ.NN)AB=-AB
        IF (ICOLBREI .EQ. 1) THEN
           INTERACT = INTERACT1
           IF(INTERACT .GT. 0) RETURN
!           BB=A1*AB*DBLE(IFAZFRCS)
!           CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI .EQ. 2) THEN
           INTERACT = 1
           RETURN
!          N=(KRA-ND1)/2+1
!          IF(((KRA-ND1)/2)*2 .EQ. (KRA-ND1)) THEN
!            CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
!            IF(DABS(S(1)).GT.EPS) THEN
!              BB=-S(1)*AB
!              IF(DABS(BB).GT.EPS)                                &
!              CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
!            END IF
!          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EL31INT
