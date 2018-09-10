!*******************************************************************
!                                                                  *
      SUBROUTINE EL41INT(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD,  &
                                                 ICOLBREI,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 10  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!     WHEN IREZ = 1   . . . . . . . . . . . . . . . . . . .        *
!                                              N'1 = N1 - 1        *
!                                              N'2 = N2 - 1        *
!                                              N'3 = N3 + 2,       *
!     WHEN IREZ = 2   . . . . . . . . . . . . . . . . . . .        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,EILE,ITREXG,IXJTIK,                *
!                        JFAZE,PERKO2,RECO,REC3                    *
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
      USE eile_I
      USE reco_I
      USE rec3_I
      USE perko2_I
      USE snrc_I
!      USE jfaze_I
      USE ixjtik_I
      USE itrexg_I
      USE itrig_I
      USE coulom_I
!      USE gg1233_I
!      USE sixj_I
!      USE speak_I
!      USE cxk_I
!      USE talk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB, &
                             JJC,JJD,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!      DIMENSION J(3)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IA,IB,IC,II,IIA,IIB,IIC,IBRD,IBRE,IP1,IG1,IP2,IG2, &
                  IAT,IID,IKK,IFAZ,IFAZP,IFAZFRCS,INN,I2,I3,JB1,JAA, &
                  JBB,JCC,J12,KRA,L1,L2,L3,ND1,ND2,NE1,NE2,N,NN,NU,  &
                  NUP1,MU
      INTEGER :: INTERACT1, INTERACT2
      INTEGER, DIMENSION(3) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
!      REAL(DOUBLE)          :: AA,AB,A1,BB,QM1,QM2,QM3,QM4,SI,RECC
      REAL(DOUBLE)          :: AA,A1,BB,SI,RECC
!      REAL(DOUBLE), DIMENSION(12) :: S
!      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      INTERACT=0
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      IF(NPEEL.LE.1)RETURN
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
!      IF(IREZ.EQ.1) THEN
!        QM1=-HALF
!        QM2=-HALF
!        QM3=HALF
!        QM4=HALF
!      ELSE
!        QM1=HALF
!        QM2=HALF
!        QM3=-HALF
!        QM4=-HALF
!      END IF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
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
!        DO II=1,20
!          COND(1,II) =ZERO
!          COND(2,II) =ZERO
!          COND(3,II) =ZERO
!          COND(4,II) =ZERO
!          COND(5,II) =ZERO
!          COND(6,II) =ZERO
!          COND(7,II) =ZERO
!          COND(8,II) =ZERO
!          COND(9,II) =ZERO
!          COND(10,II)=ZERO
!          COND(11,II)=ZERO
!          COND(12,II)=ZERO
!        END DO
      END IF
!      IFAZP=JFAZE(JC,JA,JB,JC)
!      IFAZFRCS = 1
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
!      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)+ &
!      IK3(5)*IK3(4)-ID3(5)*ID3(4)
!      IF((IFAZ/4)*4.NE.IFAZ)IFAZFRCS=-IFAZFRCS
!
!      IF(JA.GT.JB) THEN
!        JAA=JB
!        JBB=JA
!      ELSE
!        JAA=JA
!        JBB=JB
!      END IF
!      NN=0
!      JB1=JBB-1
!      DO II=JAA,JB1
!        INN=JLIST(II)
!        NN=NQ1(INN)+NN
!      END DO
!      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
! * * *                      * * *                      * * *
!     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
!           3321                                1233
!                                                    (IREZ = 1)
!     OR
!     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
!           2133                                1233
!                                                    (IREZ = 2)
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
!
        IF (ICOLBREI .EQ. 1) THEN
          INTERACT1 = 0
          IF(IREZ.EQ.2) THEN
            CALL COULOM(L1,L2,L3,L3,ID1(5),ID2(5),ID3(5),ID3(5),KRA,A1)
          ELSE
            CALL COULOM(L3,L3,L1,L2,ID3(5),ID3(5),ID1(5),ID2(5),KRA,A1)
          END IF
          IF(DABS(A1).LT.EPS) CYCLE
          INTERACT1 = 1
        END IF
!
!        AB=ZERO
        INTERACT2 = 0
        DO I3=IP1,IG1,2
          J12=(I3-1)/2
          IFAZ=J(2)-J12+1
          IF(IREZ.EQ.2)IFAZ=J(1)-J12+1
          IF((IFAZ/2)*2.NE.IFAZ) CYCLE
          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,0,IAT,AA)
          IF(IAT.EQ.0) CYCLE
          IF(IXJTIK(J(3),J(1),KRA*2,J(2),J(3),J12*2).EQ.0) CYCLE
          INTERACT2 = INTERACT2 +1
!          CALL GG1233(IK1,IK2,IK3,BK1,BK2,BK3,ID1,ID2,ID3,BD1,    &
!                    BD2,BD3,J12,QM1,QM2,QM3,QM4,AA)
!          IF(DABS(AA).LT.EPS) CYCLE
!          CALL REC3(JA,JB,JC,J(1),J(2),J12*2,1,IAT,RECC)
!          AA=AA*RECC
!          CALL SIXJ(J(3),J(1),KRA*2,J(2),J(3),J12*2,0,SI)
!          AA=AA*SI*DSQRT(DBLE(I3))
!          IFAZ=J(3)+J(1)+2*J12+2*KRA
!          IF(IREZ.EQ.2)IFAZ=J(2)+J(3)+2*J12+2*KRA
!          IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
!          AB=AB+AA
        END DO
!CIA GG KLAIDA !!!
        IF(INTERACT2 .EQ. 0) CYCLE
!        IF(DABS(AB).LT.EPS) CYCLE
!        AB=-AB*DBLE(IFAZP)
        IF (ICOLBREI .EQ. 1) THEN
           INTERACT = INTERACT1
           IF(INTERACT .GT. 0) RETURN
!          BB=A1*AB*DBLE(IFAZFRCS)
!          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI .EQ. 2) THEN
           INTERACT = 1
           RETURN
!          NU=KRA
!          IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
!            IF((ITRIG(KS(1),KS(3),NU+NU+1).NE.0) .AND.  &
!               (ITRIG(KS(2),KS(4),NU+NU+1).NE.0)) THEN
!              N=(NU-ND1)/2+1
!              IF(NU .GT. 0) THEN
!                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
!                DO MU = 1,4
!                  COND(MU,N)=COND(MU,N)+AB*S(MU)
!                END DO
!              END IF
!            END IF
!          END IF
!          NU=KRA+1
!          IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
!            IF((ITRIG(KS(1),KS(3),NU+NU-1).NE.0) .AND.  &
!               (ITRIG(KS(2),KS(4),NU+NU-1).NE.0)) THEN
!              N=(NU-ND1)/2+1
!              IF(N .LE. ND2) THEN
!                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
!                DO MU = 1,4
!                  COND(MU,N)=COND(MU,N)+AB*S(MU)
!                END DO
!              END IF
!            END IF
!          END IF
!          NU=KRA-1
!          IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
!            IF((ITRIG(KS(1),KS(3),NU+NU+3).NE.0) .AND.  &
!               (ITRIG(KS(2),KS(4),NU+NU+3).NE.0)) THEN
!              IF(NU .GE. 0) THEN
!                N=(NU-ND1)/2+1
!                IF(N .LT. ND2) THEN
!                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
!                  DO MU = 1,12
!                    COND(MU,N)=COND(MU,N)+AB*S(MU)
!                  END DO
!                END IF
!              END IF
!            END IF
!          END IF
        END IF
      END DO
!      IF (ICOLBREI .EQ. 2) THEN
!        DO N = 1,ND2
!          NU=ND1+2*(N-1)
!          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
!          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
!          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
!          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
!          IF(N.EQ.ND2) CYCLE
!          NUP1=NU+1
!          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
!          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
!        END DO
!      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL41INT  PMGG RAGG')
      STOP
      END SUBROUTINE EL41INT
