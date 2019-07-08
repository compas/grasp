!***********************************************************************
!                                                                      *
      SUBROUTINE START(IORB, ITYPE, P0, P, Q0, Q)
!                                                                      *
!   This subroutine sets up  P(1:6), Q(1:6),  required  to start the   *
!   integration for programs  OUT  and  SBSTEP .                       *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      IORB : (Input) Index of the orbital                             *
!      ITYPE: (Input) 1 = homogeneous equation; 2 = inhomogeneous      *
!             equation; 3 = variational equation                       *
!      P0   : (Input) Slope parameter                                  *
!      P    : (Output) P(1:6) are tabulated by this program            *
!      Q0   : (Output) First term in the series expansion of Q         *
!      Q    : (Output) Q(1:6) are tabulated by this program            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:48   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE CNC_C
      USE DEF_C,           ONLY: C, Z, ACCY
      USE GRID_C
      USE NPAR_C
      USE ORB_C
      USE POTE_C,          ONLY: YP, XP, XQ
      USE SCF_C,           ONLY: NDCOF, NDA, DA
      USE WAVE_C,          ONLY: PZ,PF,QF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IORB
      INTEGER, INTENT(IN) :: ITYPE
      REAL(DOUBLE), INTENT(IN) :: P0
      REAL(DOUBLE), INTENT(OUT) :: Q0
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(INOUT) :: P
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(INOUT) :: Q
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MXITER=36, KIORB, I, JORB, NITER, J
      REAL(DOUBLE), DIMENSION(6) :: RDP, RDQ, RSEP, RSEQ
      REAL(DOUBLE), DIMENSION(2:6) :: SPEST, SQEST
      REAL(DOUBLE) :: OBC, ZONC, GIORB, FKIORB, OMGI, OMGMK, OMGPK, RSEP1, &
         RSEQ1, PZERO, P1, Q1, SUMP, SUMQ, FACTOR, CSQ, TWOCSQ, ENERGY, ENEFAC&
         , RI, RPI, RIRPI, YPIRPI, DIFMAW, DIFMAX, COEFIJ, PI, QI, RJ
!-----------------------------------------------
!
!
!   Initialization
!
      OBC = 1.0D00/C
      ZONC = Z*OBC
      GIORB = GAMA(IORB)
      KIORB = NAK(IORB)
      FKIORB = DBLE(KIORB)
      OMGI = 1.0D00 - GIORB
!
      OMGMK = OMGI - FKIORB
      OMGPK = OMGI + FKIORB
!
!   Determine P(1), Q(1): THESE STORE  R**(-GAMMA)*(P(1),Q(1));
!   set up  RSEP  and  RSEQ , the inhomogeneous terms
!
      IF (ITYPE==1 .OR. ITYPE==2) THEN
         IF (NPARM == 0) THEN
            P(1) = P0
            IF (KIORB < 0) THEN
               Q(1) = -P0*ZONC/(GIORB - FKIORB)
            ELSE
               Q(1) = P0*(GIORB + FKIORB)/ZONC
            ENDIF
         ELSE
            IF (KIORB < 0) THEN
               P(1) = P0
               Q(1) = 0.0D00
            ELSE
               P(1) = 0.0D00
               Q(1) = P0*(GIORB + FKIORB)/ZONC
            ENDIF
         ENDIF
         IF (ITYPE == 1) THEN
            RSEP = 0.0D00
            RSEQ = 0.0D00
         ELSE
            RSEP1 = 0.0D00
            RSEQ1 = 0.0D00
            DO I = 1, NDCOF
               JORB = NDA(I)
               PZERO = PZ(JORB)
               IF (NPARM == 0) THEN
                  P1 = PZERO
                  IF (KIORB < 0) THEN
                     Q1 = -PZERO*ZONC/(GIORB - FKIORB)
                  ELSE
                     Q1 = PZERO*(GIORB + FKIORB)/ZONC
                  ENDIF
                  SUMP = ZONC*Q1
                  SUMQ = -ZONC*P1
               ELSE
                  IF (KIORB < 0) THEN
                     P1 = PZERO
                     Q1 = 0.0D00
                  ELSE
                     P1 = 0.0D00
                     Q1 = PZERO*(GIORB + FKIORB)/ZONC
                  ENDIF
                  SUMP = 0.0D00
                  SUMQ = 0.0D00
               ENDIF
               FACTOR = DA(I)
               RSEP1 = RSEP1 + FACTOR*(SUMP + OMGMK*P1)
               RSEQ1 = RSEQ1 + FACTOR*(SUMQ + OMGPK*Q1)
            END DO
            FACTOR = RP(1)
            RSEP(1) = FACTOR*RSEP1
            RSEQ(1) = FACTOR*RSEQ1
            DO I = 2, 6
               FACTOR = -RP(I)*R(I)**(-GIORB)
               RSEP(I) = FACTOR*XP(I)
               RSEQ(I) = FACTOR*XQ(I)
            END DO
         ENDIF
      ELSE IF (ITYPE == 3) THEN
         P(1) = 0.0D00
         Q(1) = 0.0D00
         RSEP(1) = 0.0D00
         RSEQ(1) = 0.0D00
         DO I = 2, 6
            FACTOR = OBC*RP(I)*R(I)**OMGI
            RSEP(I) = -FACTOR*QF(I,IORB)
            RSEQ(I) = FACTOR*PF(I,IORB)
         END DO
      ENDIF
      Q0 = Q(1)
!
!   Set up  RDP  and  RDQ
!
      CSQ = C*C
      TWOCSQ = CSQ + CSQ
      ENERGY = E(IORB)
      ENEFAC = TWOCSQ - ENERGY
      DO I = 1, 6
         RI = R(I)
         RPI = RP(I)
         RIRPI = RI*RPI
         YPIRPI = YP(I)*RPI
         RDP(I) = -OBC*(ENEFAC*RIRPI + YPIRPI)
         RDQ(I) = -OBC*(ENERGY*RIRPI - YPIRPI)
      END DO
!
!   Determine  P(2:6) , Q(2:6)
!
!   Initilizations for the iterations
!
      NITER = 0
      P1 = P(1)
      Q1 = Q(1)
      DIFMAW = MAX(ABS(P1),ABS(Q1))
!
      P(2:6) = P1
      Q(2:6) = Q1
!
!   This is the largest factor by which any result will be
!   multiplied
!
      FACTOR = R(6)**GIORB
!
!   Now iterate
!
    7 CONTINUE
      NITER = NITER + 1
      DIFMAX = 0.0D00
      DO J = 2, 6
         SUMP = SUM(CNC6C(:,J)*(OMGMK*RP(:6)*P(:6)-RDP*Q(:6)+RSEP))
         SUMQ = SUM(CNC6C(:,J)*(OMGPK*RP(:6)*Q(:6)-RDQ*P(:6)+RSEQ))
         RJ = R(J)
         SUMP = SUMP/RJ
         SUMQ = SUMQ/RJ
         SPEST(J) = SUMP
         SQEST(J) = SUMQ
         DIFMAX = MAX(DIFMAX,ABS(SUMP - P(J)))
         DIFMAX = MAX(DIFMAX,ABS(SUMQ - Q(J)))
      END DO
!zou  IF (DIFMAX .LT. DIFMAW) THEN
      P(2:6) = SPEST(:6)
      Q(2:6) = SQEST(:6)
      DIFMAW = DIFMAX
      DIFMAX = DIFMAX*FACTOR
      IF (DIFMAX > ACCY) THEN
         IF (NITER < MXITER) THEN
            GO TO 7
         ENDIF
      ENDIF
!     ELSE
!        DIFMAX = DIFMAX*FACTOR
!        IF (DIFMAX .GT. ACCY) THEN
!           WRITE (*,300) NP(IORB),NH(IORB),DIFMAX,NITER,ACCY
!        ENDIF
!zou  ENDIF
! not convergent, using the initial P,Q
      IF (DIFMAX > ACCY) THEN
         P(2:6) = P1
         Q(2:6) = Q1
      ENDIF
!zou
!
!   All done
!
!   This is always true in GRASP2
!
      P(1) = 0.0D00
      Q(1) = 0.0D00
!
      DO I = 2, 6
         FACTOR = R(I)**GIORB
         P(I) = FACTOR*P(I)
         Q(I) = FACTOR*Q(I)
      END DO
!
      RETURN
!
  300 FORMAT('START: ',1I2,1A2,' subshell: accuracy ',1P,1D8.1,/,&
         ' attained after ',1I2,' iterations; this fails the'/,&
         ' accuracy criterion ',D8.1,'.')
      RETURN
!
      END SUBROUTINE START
