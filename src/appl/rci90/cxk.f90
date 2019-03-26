!***********************************************************************
!                                                                      *
      SUBROUTINE CXK(S, IS, KAPS, NU, K, IBR, IEX)
!                                                                      *
!   Computes  the  coefficients of radial integrals in the expansion   *
!   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS  : Orbital labels                                            *
!      KAPS: Values of 2*kappa                                         *
!      NU  : Order of radial integral                                  *
!      K   : Index of tensor operator                                  *
!      IEX : 1 for direct, 2 for exchange terms                        *
!      IBR : Classifies type of radial integral.There are 4 distinct   *
!            cases:                                                    *
!            IBR = 1 A. All states distinct                            *
!                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
!                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
!                    These give 12 distinct  radial  integrals, with   *
!                    values of K and NU limited only by angular mom-   *
!                    entum and parity                                  *
!            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
!                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
!                    This case gives one non-zero integral when K =    *
!                    NU is ODD                                         *
!            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
!                    (IA .NE. IB)                                      *
!                    Integrals of magnetic F-type when K = NU is odd   *
!            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
!                    netic G-type integrals and  four  H-TYPE  inte-   *
!                    grals                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!      S   : Coefficients S(MU) MU = 1,12                              *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92] CRE.                                           *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cre_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NU
      INTEGER  :: K
      INTEGER, INTENT(IN) :: IBR
      INTEGER, INTENT(IN) :: IEX
      INTEGER, INTENT(IN) :: IS(4)
      INTEGER, INTENT(IN) :: KAPS(4)
      REAL(DOUBLE), INTENT(INOUT) :: S(12)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MU, IA, IB, IC, ID, KA, KB, KC, KD, KK, IK, IP
      REAL(DOUBLE) :: D, H, DK1, DK2, FK, GK, G1, G2, G3, G4, A, F1, F2, F3, F4&
         , B, DK
!-----------------------------------------------
!
!
!   1.0  Initialization
!
      S = 0.0D00
!
      IA = IS(1)
      IB = IS(2)
      IC = IS(3)
      ID = IS(4)
      KA = KAPS(1)/2
      KB = KAPS(2)/2
      KC = KAPS(3)/2
      KD = KAPS(4)/2
      IF (IEX == 2) THEN
         KK = KD
         IK = ID
         KD = KC
         ID = IC
         KC = KK
         IC = IK
      ENDIF
      SELECT CASE (IBR)
      CASE DEFAULT
         GO TO 17
!
!   2.0  IBR = 1 --- The general case
!
      CASE (1)
         IF (NU - K >= 0) THEN
            IF (NU - K <= 0) THEN
               S(1) = -(KA + KC)*(KD + KB)
               IF (K == 0) GO TO 16
               D = K*(K + 1)
               H = CRE(KA,K,KC)*CRE(KB,K,KD)
               IF (MOD(K,2) /= 0) H = -H
               S(1) = S(1)*H/D
               S(2:4) = S(1)
               RETURN
            ENDIF
!
!   2.2  NU = K+1
!
            DK1 = KC - KA
            DK2 = KD - KB
            FK = K
            GK = K + 1
            G1 = DK1 - GK
            G2 = DK1 + GK
            G3 = DK2 - GK
            G4 = DK2 + GK
            KK = K + K + 1
            H = CRE(KA,K,KC)*CRE(KB,K,KD)
            IF (MOD(K,2) /= 0) H = -H
            A = H*FK/GK/DBLE(KK*(KK + 2))
            S(1) = A*G1*G3
            S(2) = A*G2*G4
            S(3) = A*G1*G4
            S(4) = A*G2*G3
            RETURN
         ENDIF
!
!   2.2  NU = K-1
!
         DK1 = KC - KA
         DK2 = KD - KB
         FK = K
         GK = K + 1
         F1 = DK1 - FK
         F2 = DK1 + FK
         F3 = DK2 - FK
         F4 = DK2 + FK
         G1 = DK1 - GK
         G2 = DK1 + GK
         G3 = DK2 - GK
         G4 = DK2 + GK
         KK = K + K + 1
         H = CRE(KA,K,KC)*CRE(KB,K,KD)
         IF (MOD(K,2) /= 0) H = -H
         A = H*GK/FK/DBLE(KK*(KK - 2))
         S(1) = A*F2*F4
         S(2) = A*F1*F3
         S(3) = A*F2*F3
         S(4) = A*F1*F4
         B = H/DBLE(KK*KK)
         S(5) = B*F2*G3
         S(6) = B*F4*G1
         S(7) = B*F1*G4
         S(8) = B*F3*G2
         S(9) = B*F2*G4
         S(10) = B*F3*G1
         S(11) = B*F1*G3
         S(12) = B*F4*G2
         RETURN
!
!   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
!
      CASE (2)
         IF (IA/=IC .OR. IB==ID) THEN
            IF (IA==IC .OR. IB/=ID) GO TO 17
!
            IK = IB
            IB = IA
            IA = IK
            IK = ID
            ID = IC
            IC = IK
!
            KK = KB
            KB = KA
            KA = KK
            KK = KD
            KD = KC
            KC = KK
         ENDIF
!
         IF (MOD(K,2) /= 1) RETURN
         DK = K*(K + 1)
         H = CRE(KA,K,KC)*CRE(KB,K,KD)/DK
         S(1) = H*DBLE(4*KA*(KB + KD))
         RETURN
!
!   4.0  IBR = 3. Direct magnetic F-integrals
!
      CASE (3)
         IF (IA/=IC .OR. IB/=ID) GO TO 17
         IF (MOD(K,2) /= 1) RETURN
         DK = K*(K + 1)
         H = CRE(KA,K,KA)*CRE(KB,K,KB)/DK
         S(1) = H*DBLE(16*KA*KB)
         RETURN
!
!   5.0   IBR = 4. Exchange magnetic G- and H-integrals
!
      CASE (4)
         IF (IA/=ID .OR. IB/=IC) GO TO 17
         IF (NU - K >= 0) THEN
            IF (NU - K <= 0) THEN
               S(1) = DBLE(KA + KB)*CRE(KA,K,KB)
               IP = ABS(KA) - ABS(KB) + K + 1
               S(1) = S(1)*S(1)/DBLE(K*(K + 1))
               IF (MOD(IP,2) /= 0) S(1) = -S(1)
               S(3) = S(1)
               S(2) = S(1) + S(1)
               RETURN
            ENDIF
!
!   5.2  NU = K+1
!
            DK = KB - KA
            GK = K + 1
            FK = K
            G1 = DK + GK
            G2 = DK - GK
            KK = K + K + 1
            H = CRE(KA,K,KB)**2
            IF (KA*KB < 0) H = -H
            A = H*FK/GK/DBLE(KK*(KK + 2))
            S(1) = -A*G1*G1
            S(2) = -2.0D00*A*G1*G2
            S(3) = -A*G2*G2
            RETURN
         ENDIF
!
!   5.3  NU = K-1
!
         DK = KB - KA
         FK = K
         GK = K + 1
         F1 = DK + FK
         F2 = DK - FK
         G1 = DK + GK
         G2 = DK - GK
         KK = K + K + 1
         H = CRE(KA,K,KB)**2
         IF (KA*KB < 0) H = -H
         A = H*GK/FK/DBLE(KK*(KK - 2))
         S(1) = -A*F2*F2
         S(2) = -2.0D00*A*F1*F2
         S(3) = -A*F1*F1
         B = H/DBLE(KK*KK)
         B = B + B
         S(4) = -B*F1*G2
!     S(5) = S(4)
         S(5) = -B*F2*G1
         S(6) = -B*F1*G1
         S(7) = -B*F2*G2
         RETURN
      END SELECT
!
!   6.0  Special cases and errors
!
!   Illegal zero value of K in Type 1
!
   16 CONTINUE
      WRITE (*, 300) IS(1), IS(2), IS(3), IS(4), NU, IBR, IEX
      STOP
!
!   Illegal combination of states in Type 3 or 4
!
   17 CONTINUE
      WRITE (*, 301) IBR, IS(1), IS(2), IS(3), IS(4), NU, K, IEX
      STOP
!
  300 FORMAT('CXK: Illegal value K = 0 -'/,1X,4I3,2X,I3,2X,2I2)
  301 FORMAT('CXK: Type ',I2,'-'/,1X,I2,3X,4I3,2X,2I3,2X,I2)
      RETURN
!
      END SUBROUTINE CXK
