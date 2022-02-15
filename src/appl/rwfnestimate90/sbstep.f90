

!***********************************************************************
!                                                                      *
      SUBROUTINE SBSTEP(IORB, NSTRT, NEND, P, Q)
!                                                                      *
!   This  subroutine continues the solution of the homogeneous Dirac   *
!   radial equation  from tabulation point NSTRT to tabulation point   *
!   NEND. The algorithm of J E Sienkiewicz and W E Baylis, J Phys B:   *
!   At Mol Phys 20 (1987) 5145-5156, p 5155, is used.                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def, ONLY: NNNP, NNN1, NNNW
      !USE DEF_C
      USE GRID_C
      USE INT_C, ONLY: TF, TG
      USE ORB_C
      USE SBC_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IORB
      INTEGER , INTENT(IN) :: NSTRT
      INTEGER , INTENT(IN) :: NEND
      REAL(DOUBLE) , INTENT(INOUT) :: P(NNNP)
      REAL(DOUBLE) , INTENT(INOUT) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: B1 = -0.50D00
      REAL(DOUBLE), PARAMETER :: B2 = 0.50D00
      REAL(DOUBLE), PARAMETER :: B3 = 0.75D00
      REAL(DOUBLE), PARAMETER :: B4 = 0.25D00
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IDIFF, LOC, J, JP1, JM1
      REAL(DOUBLE) :: TBH, TC1, FK, CC, PJ, QJ, FAC, PPJ, QPJ, PJM1, QJM1, &
         PPJM1, QPJM1, PJM2, QJM2, PPJM2, QPJM2, PJM3, QJM3, PPJM3, QPJM3, PJM4&
         , QJM4, PPJM4, QPJM4, RPPJ, RQPJ, CCRPOR, CPJP1, CMJP1, FJP1, GJP1, &
         DENOM, FACTOR, PJP1, QJP1, PPJP1, QPJP1, PJP2, QJP2, PPJP2, QPJP2, &
         PJP3, QJP3, PPJP3, QPJP3, PJP4, QJP4, PPJP4, QPJP4, RPMJ, RQMJ, CPJM1&
         , CMJM1, FJM1, GJM1
      REAL(DOUBLE) :: C1, C2, C3, C4, C5, C6
      LOGICAL :: OUT
!-----------------------------------------------
!
!
!   Initializations
!
      C1 = C(1)
      C2 = C(2)
      C3 = C(3)
      C4 = C(4)
      C5 = C(5)
      C6 = C(6)
      TBH = 2.0D00/H
      TC1 = 2.0D00*C1
!
!   Determine whether integration is inward or outwarD
!
      IDIFF = NEND - NSTRT
      IF (IDIFF > 0) THEN
         OUT = .TRUE.
      ELSE IF (IDIFF == 0) THEN
         RETURN
      ELSE IF (IDIFF < 0) THEN
         OUT = .FALSE.
      ENDIF
!
!   Overall initializations
!
      FK = DBLE(NAK(IORB))
      CC = C1*FK*H
!
!   Perform integration depending on case
!
      IF (OUT) THEN
!
!   Initialization for outward integration
!
         LOC = NSTRT
         PJ = P(LOC)
         QJ = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJ = (-FAC*PJ) - TBH*TF(LOC)*QJ
         QPJ = FAC*QJ - TBH*TG(LOC)*PJ
!
         LOC = LOC - 1
         PJM1 = P(LOC)
         QJM1 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM1 = (-FAC*PJM1) - TBH*TF(LOC)*QJM1
         QPJM1 = FAC*QJM1 - TBH*TG(LOC)*PJM1
!
         LOC = LOC - 1
         PJM2 = P(LOC)
         QJM2 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM2 = (-FAC*PJM2) - TBH*TF(LOC)*QJM2
         QPJM2 = FAC*QJM2 - TBH*TG(LOC)*PJM2
!
         LOC = LOC - 1
         PJM3 = P(LOC)
         QJM3 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM3 = (-FAC*PJM3) - TBH*TF(LOC)*QJM3
         QPJM3 = FAC*QJM3 - TBH*TG(LOC)*PJM3
!
         LOC = LOC - 1
         PJM4 = P(LOC)
         QJM4 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM4 = (-FAC*PJM4) - TBH*TF(LOC)*QJM4
         QPJM4 = FAC*QJM4 - TBH*TG(LOC)*PJM4
!
!   March out
!
         J = NSTRT - 1
    1    CONTINUE
         J = J + 1
!
         RPPJ = B1*PJ + B2*PJM1 + B3*PJM2 + B4*PJM3 + C2*PPJ + C3*PPJM1 + C4*&
            PPJM2 + C5*PPJM3 + C6*PPJM4
         RQPJ = B1*QJ + B2*QJM1 + B3*QJM2 + B4*QJM3 + C2*QPJ + C3*QPJM1 + C4*&
            QPJM2 + C5*QPJM3 + C6*QPJM4
!
         JP1 = J + 1
         CCRPOR = CC*RPOR(JP1)
         CPJP1 = 1.0D00 + CCRPOR
         CMJP1 = 1.0D00 - CCRPOR
         FJP1 = TC1*TF(JP1)
         GJP1 = TC1*TG(JP1)
         DENOM = CPJP1*CMJP1 - GJP1*FJP1
         FACTOR = 1.0D00/DENOM
         PJP1 = (CMJP1*RPPJ - FJP1*RQPJ)*FACTOR
         QJP1 = (CPJP1*RQPJ - GJP1*RPPJ)*FACTOR
         P(JP1) = PJP1
         Q(JP1) = QJP1
!
         IF (JP1 < NEND) THEN
!
            PPJM4 = PPJM3
            QPJM4 = QPJM3
!
            PJM3 = PJM2
            QJM3 = QJM2
            PPJM3 = PPJM2
            QPJM3 = QPJM2
!
            PJM2 = PJM1
            QJM2 = QJM1
            PPJM2 = PPJM1
            QPJM2 = QPJM1
!
            PJM1 = PJ
            QJM1 = QJ
            PPJM1 = PPJ
            QPJM1 = QPJ
!
            PJ = PJP1
            QJ = QJP1
            FAC = FK*RPOR(JP1)
            PPJ = (-FAC*PJ) - TBH*TF(JP1)*QJ
            QPJ = FAC*QJ - TBH*TG(JP1)*PJ
!
            GO TO 1
!
         ENDIF
      ELSE
!
!   Initializations for inward integration
!
         LOC = NSTRT
         PJ = P(LOC)
         QJ = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJ = (-FAC*PJ) - TBH*TF(LOC)*QJ
         QPJ = FAC*QJ - TBH*TG(LOC)*PJ
!
         LOC = LOC + 1
         PJP1 = P(LOC)
         QJP1 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP1 = (-FAC*PJP1) - TBH*TF(LOC)*QJP1
         QPJP1 = FAC*QJP1 - TBH*TG(LOC)*PJP1
!
         LOC = LOC + 1
         PJP2 = P(LOC)
         QJP2 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP2 = (-FAC*PJP2) - TBH*TF(LOC)*QJP2
         QPJP2 = FAC*QJP2 - TBH*TG(LOC)*PJP2
!
         LOC = LOC + 1
         PJP3 = P(LOC)
         QJP3 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP3 = (-FAC*PJP3) - TBH*TF(LOC)*QJP3
         QPJP3 = FAC*QJP3 - TBH*TG(LOC)*PJP3
!
         LOC = LOC + 1
         PJP4 = P(LOC)
         QJP4 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP4 = (-FAC*PJP4) - TBH*TF(LOC)*QJP4
         QPJP4 = FAC*QJP4 - TBH*TG(LOC)*PJP4
!
!   March in
!
         J = NSTRT + 1
    2    CONTINUE
         J = J - 1
!
         RPMJ = B1*PJ + B2*PJP1 + B3*PJP2 + B4*PJP3 - C2*PPJ - C3*PPJP1 - C4*&
            PPJP2 - C5*PPJP3 - C6*PPJP4
         RQMJ = B1*QJ + B2*QJP1 + B3*QJP2 + B4*QJP3 - C2*QPJ - C3*QPJP1 - C4*&
            QPJP2 - C5*QPJP3 - C6*QPJP4
!
         JM1 = J - 1
         CCRPOR = CC*RPOR(JM1)
         CPJM1 = 1.0D00 + CCRPOR
         CMJM1 = 1.0D00 - CCRPOR
         FJM1 = TC1*TF(JM1)
         GJM1 = TC1*TG(JM1)
         DENOM = CPJM1*CMJM1 - GJM1*FJM1
         FACTOR = 1.0D00/DENOM
         PJM1 = (CPJM1*RPMJ + FJM1*RQMJ)*FACTOR
         QJM1 = (CMJM1*RQMJ + GJM1*RPMJ)*FACTOR
         P(JM1) = PJM1
         Q(JM1) = QJM1
!
         IF (JM1 > NEND) THEN
!
            PPJP4 = PPJP3
            QPJP4 = QPJP3
!
            PJP3 = PJP2
            QJP3 = QJP2
            PPJP3 = PPJP2
            QPJP3 = QPJP2
!
            PJP2 = PJP1
            QJP2 = QJP1
            PPJP2 = PPJP1
            QPJP2 = QPJP1
!
            PJP1 = PJ
            QJP1 = QJ
            PPJP1 = PPJ
            QPJP1 = QPJ
!
            PJ = PJM1
            QJ = QJM1
            FAC = FK*RPOR(JM1)
            PPJ = (-FAC*PJ) - TBH*TF(JM1)*QJ
            QPJ = FAC*QJ - TBH*TG(JM1)*PJ
!
            GO TO 2
!
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE SBSTEP
