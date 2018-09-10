!********************************************************************
!                                                                  *
      SUBROUTINE EL2INT(JJA,JJB,JA,JB,ICOLBREI,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
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
      USE itrexg_I
      USE itrig_I
      USE ixjtik_I
      USE snrc_I
      USE coulom_I
!      USE gg1122_I
!      USE sixj_I
!      USE cxk_I
!      USE talk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IAT,IA,IB,IBRD,IBRE,IP1,IP2,IG1,IG2,IKK,II,I2,I3, &
!                 IFAZ,IFAZ1,IFAZFRCS,J12,JAA,JBB, &
                 IFAZ,IFAZ1,J12,JAA,JBB, &
                 KRA,L1,L2,N,ND1,ND2,NE1,NE2,NUP1,NU,MU
      INTEGER :: INTERACT1, INTERACT2
      INTEGER, DIMENSION(2) :: J
      INTEGER, DIMENSION(4) :: IS,KAPS,KS
!      REAL(DOUBLE)          :: AA,A1,AB,BB,QM1,QM2,QM3,QM4,RECC,SI
      REAL(DOUBLE)          :: AA,A1,BB,RECC
!      REAL(DOUBLE), DIMENSION(12)    :: S
!      REAL(DOUBLE), DIMENSION(12,20) :: COND
!-----------------------------------------------
      INTERACT = 0
      IF(NPEEL.LE.1)RETURN
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
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IB
        IS(4)=IB
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
      DO I2=IP2,IG2,2
        KRA=(I2-1)/2
        IF (ICOLBREI .EQ. 1) THEN
          INTERACT1 = 0
          CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
          IF(DABS(A1).LT.EPS) CYCLE
          INTERACT1 = 1
!          A1=-HALF*A1
        END IF
!        AB=ZERO
          INTERACT2 = 0
          DO I3=IP1,IG1,2
            J12=(I3-1)/2
            CALL RECO2(JAA,JBB,J12*2,0,IAT,RECC)
            IF(IAT.NE.0) THEN
              IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2).NE.0) THEN
                 INTERACT2 = INTERACT2 + 1
!                CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
!                IF(DABS(AA).GT.EPS) THEN
!                  CALL RECO2(JAA,JBB,J12*2,1,IAT,RECC)
!                  AA=AA*RECC
!                  CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
!                  AA=AA*SI*DSQRT(DBLE(I3))
!                  IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
!                  IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
!                  AB=AB+AA
!                END IF
              END IF
            END IF
          END DO
          IF(INTERACT2 .EQ. 0) CYCLE
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
!        IFAZFRCS=1
!        IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
!        IF((IFAZ1/4)*4.NE.IFAZ1)IFAZFRCS=-IFAZFRCS
!
        IF (ICOLBREI .EQ. 1) THEN
           INTERACT = INTERACT1
           IF(INTERACT .GT. 0) RETURN
!          BB=A1*AB*DBLE(IFAZFRCS)
!          IF(DABS(BB).GT.EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
        ELSE IF (ICOLBREI .EQ. 2) THEN
           INTERACT = 1
           RETURN
!          NU=KRA
!          IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
!            IF((ITRIG(KS(1),KS(3),NU+NU+1).NE.0) .AND.   &
!               (ITRIG(KS(2),KS(4),NU+NU+1).NE.0)) THEN
!              N=(NU-ND1)/2+1
!              IF(NU .GT. 0) THEN
!                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
!                DO MU = 1,4
!                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
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
!                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
!                END DO
!              END IF
!            END IF
!          END IF
!          NU=KRA-1
!          IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
!            IF((ITRIG(KS(1),KS(3),NU+NU+3).NE.0) .AND.   &
!               (ITRIG(KS(2),KS(4),NU+NU+3).NE.0)) THEN
!              IF(NU .GE. 0) THEN
!                N=(NU-ND1)/2+1
!                IF(N .LT. ND2) THEN
!                  CALL CXK(S,IS,KAPS,NU,KRA,1,1)
!                  DO MU = 1,12
!                    COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
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
!          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
!          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
!          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
!          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
!          IF(N.EQ.ND2) CYCLE
!          NUP1=NU+1
!          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
!          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
!          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
!          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
!          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
!          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
!          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
!          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
!        END DO
!      END IF
      RETURN
      END SUBROUTINE EL2INT
