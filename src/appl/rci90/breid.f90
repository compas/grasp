!***********************************************************************
!                                                                      *
      SUBROUTINE BREID(JA, JB, JA1, IPCA, JB1)
!                                                                      *
!   Computes closed shell contributions - aaaa and exchange only.      *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
!                                                                      *
!                                           LAST UPDATE: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE bcore_C
      USE cons_C
      USE debug_C
      USE m_C
      USE orb_C, ONLY: np, nak
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE snrc_I
      USE clrx_I
      USE talk_I
      USE itrig_I
      USE cxk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JA
      INTEGER  :: JB
      INTEGER, INTENT(IN) :: JA1
      INTEGER, INTENT(IN) :: IPCA
      INTEGER, INTENT(IN) :: JB1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NUMAX = 20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: JS, KAPS, KS
      INTEGER :: IA1, IB1, ISG, NQS1, NQS2, I, ND1, ND2, NE1, NE2, IBRD, IBRE, &
         N, NU, K, KAP1, ITYPE, MU, IP, IPP, KK, NUP1
      REAL(DOUBLE), DIMENSION(7,20) :: CONE
      REAL(DOUBLE), DIMENSION(12) :: S
      REAL(DOUBLE) :: CONST, GAM, DKSKS, DNUNU1, COEF, PROC, PROD
!-----------------------------------------------
!
!   1.0  Initialization
!
      IF (IPCA == 2) THEN
         IA1 = KLIST(JA1)
      ELSE
         IA1 = JLIST(JA1)
      ENDIF
      IB1 = KLIST(JB1)
!
      ISG = 1
      IF (JA == JB) THEN
         IF (ICORE(IA1)/=0 .AND. ICORE(IB1)/=0) THEN
            IF (JA > 1) RETURN
            ISG = -1
         ENDIF
      ENDIF
!
      JS(1) = IA1
      JS(2) = IB1
      JS(3) = IA1
      JS(4) = IB1
      NQS1 = NQ1(IA1)
      NQS2 = NQ2(IB1)
      DO I = 1, 4
         KAPS(I) = 2*NAK(JS(I))
         KS(I) = IABS(KAPS(I))
      END DO
      CONST = NQS1*NQS2
      IF (IBUG2 /= 0) WRITE (99, 300) IA1, IB1
!
!   2.0  Set range of tensor indices
!
      CALL SNRC (JS, KAPS, KS, ND1, ND2, NE1, NE2, IBRD, IBRE)
      IF (IBUG2 /= 0) WRITE (99, 301) ND1, ND2, NE1, NE2, IBRD, IBRE
      IF (IA1 == IB1) THEN
!
!   3.0 Calculate aaaa interaction
!
         DO N = 1, ND2
            NU = ND1 + 2*(N - 1)
            K = NU
            IF (MOD(K,2) /= 1) RETURN
            KAP1 = KAPS(1)/2
            GAM = CLRX(KAP1,NU,KAP1)
            DKSKS = KS(1)*KS(1)
            DNUNU1 = NU*(NU + 1)
            COEF = CONST*TWO*DKSKS*GAM*GAM/DNUNU1
            IF (IBUG2 /= 0) WRITE (99, 302) NU, GAM, COEF
            ITYPE = ISG*4
            CALL TALK (JA, JB, NU, IA1, IA1, IA1, IA1, ITYPE, COEF)
         END DO
         RETURN
      ENDIF
!
!   Calculate exchange interactions
!
      IF (IBRE < 0) RETURN
      IF (NE2 > NUMAX) THEN
         WRITE (*, 304)
         STOP
      ENDIF
!
      CONE(:,:NE2) = ZERO
!
      PROC = -CONST/DBLE(KS(1)*KS(2))
!
!   Negative sign arises from Pauli phase factor
!
      DO N = 1, NE2
         NU = NE1 + 2*(N - 1)
         K = NU
         IP = (KS(1)-KS(2))/2 + K
         IPP = IP + 1
         IF (NU /= 0) THEN
            KK = K + K + 1
            IF (ITRIG(KS(1),KS(2),KK) /= 0) THEN
               PROD = PROC
               IF (MOD(IP,2) /= 0) PROD = -PROD
               CALL CXK (S, JS, KAPS, NU, K, IBRE, 2)
               IF (IBUG2 /= 0) WRITE (99, 303) PROD, (S(MU),MU=1,3)
               CONE(:3,N) = CONE(:3,N) + PROD*S(:3)
            ENDIF
!
            K = NU - 1
            KK = K + K + 1
            IF (ITRIG(KS(1),KS(2),KK) /= 0) THEN
               PROD = PROC
               IF (MOD(IPP,2) /= 0) PROD = -PROD
               CALL CXK (S, JS, KAPS, NU, K, IBRE, 2)
               IF (IBUG2 /= 0) WRITE (99, 303) PROD, (S(MU),MU=1,3)
               CONE(:3,N) = CONE(:3,N) + PROD*S(:3)
!
            ENDIF
         ENDIF
         IF (N == NE2) EXIT
         K = NU + 1
         KK = K + K + 1
         PROD = PROC
         IF (MOD(IPP,2) /= 0) PROD = -PROD
         CALL CXK (S, JS, KAPS, NU, K, IBRE, 2)
         IF (IBUG2 /= 0) WRITE (99, 303) PROD, (S(MU),MU=1,7)
         CONE(:,N) = CONE(:,N) + PROD*S(:7)
      END DO
!
!   4.0  Output results
!
      DO N = 1, NE2
         NU = NE1 + 2*(N - 1)
         ITYPE = ISG*5
         CALL TALK (JA, JB, NU, IB1, IA1, IB1, IA1, ITYPE, CONE(1,N))
         CALL TALK (JA, JB, NU, IA1, IB1, IB1, IA1, ITYPE, CONE(2,N))
         CALL TALK (JA, JB, NU, IA1, IB1, IA1, IB1, ITYPE, CONE(3,N))
         IF (N == NE2) CYCLE
         NUP1 = NU + 1
         ITYPE = ISG*6
         CALL TALK (JA, JB, NUP1, IA1, IB1, IA1, IB1, ITYPE, CONE(4,N))
         CALL TALK (JA, JB, NUP1, IB1, IA1, IB1, IA1, ITYPE, CONE(5,N))
         CALL TALK (JA, JB, NUP1, IA1, IB1, IB1, IA1, ITYPE, CONE(6,N))
         CALL TALK (JA, JB, NUP1, IB1, IA1, IA1, IB1, ITYPE, CONE(7,N))
      END DO
      RETURN
!
  300 FORMAT('BREID: orbitals ',2I3)
  301 FORMAT(2X,'ND1 ND2 NE1 NE2 IBRD IBRE ',6I5)
  302 FORMAT(2X,'aaaa contribution: NU,GAM,COEF',I5,2(3X,1P,D15.8))
  303 FORMAT(2X,'PROD = ',1P,D15.8,/,' S',7D15.8)
  304 FORMAT('BREID: Dimension error for NUMAX.')
      RETURN
!
      END SUBROUTINE BREID
