

!***********************************************************************
!                                                                      *
      SUBROUTINE SOLVH(IORB, FAIL)
!                                                                      *
!   This routine solves the homogeneous Dirac radial equation.         *
!                                                                      *
!   Arguments:  IORB : (Input) Index of orbital                        *
!               FAIL : (Output) .TRUE. if solution not obtained        *
!                                                                      *
!   The direct potential is assumed tabulated in the COMMON array YP   *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RSCF92]: COUNT, SBSTEP, SETPOT, START, TAIL.          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEF_C
      USE GRID_C
      USE ORB_C
      USE POTE_C, ONLY: YP
      USE TATB_C, ONLY: TA, MTP
      USE WAVE_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setpot_I
      USE start_I
      USE sbstep_I
      USE tail_I
      USE count_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IORB
      LOGICAL , INTENT(OUT) :: FAIL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MXK, NPIORB, NKIORB, NAKABS, NLIORB, NNP, NREL, KOUNT, JP, &
         ITYPE, NSTRT, I, NPC
      REAL(DOUBLE) :: EPSLON, CSQ, ALPHA, FKABS, FKAP2, ZALPHA, GAMMA, EBYM, &
         EMIN, EMAX, DELE, EEST, Q0, PFJPO, QFJPO, PFJPI, QFJPI, RATIO, SGN, &
         DNORM, DMSMCH, DNFAC
!-----------------------------------------------
!
      DATA MXK/ 75/
!
!   A solution is deemed continuous when the relative mismatch in the
!   small component is within EPSLON
!
      EPSLON = ACCY*0.1E00
!
!   Establish the number of nodes in the large component
!
      NPIORB = NP(IORB)
      NKIORB = NAK(IORB)
      NAKABS = ABS(NKIORB)
      IF (NKIORB < 0) THEN
         NLIORB = NAKABS - 1
      ELSE
         NLIORB = NAKABS
      ENDIF
      NNP = NPIORB - NLIORB - 1
!
!   Establish the bounds on, and an estimate of, the eigenvalue
!
      CSQ = C*C
      ALPHA = 1.0D00/C
      NREL = NPIORB - NAKABS
      FKABS = DBLE(NAKABS)
      FKAP2 = FKABS*FKABS
!
      ZALPHA = YP(N)*ALPHA
      IF (ZALPHA < FKABS) THEN
         GAMMA = SQRT(FKAP2 - ZALPHA*ZALPHA)
         EBYM = 1.0D00/SQRT(1.0D00 + (ZALPHA/(GAMMA + NREL + 0.5D00))**2)
         EMIN = (1.0D00 - EBYM)*CSQ
      ELSE
         EMIN = 0.25D00*CSQ/DBLE(NPIORB*NPIORB)
      ENDIF
!
      ZALPHA = Z*ALPHA
!
      IF (ZALPHA < FKABS) THEN
         GAMMA = SQRT(FKAP2 - ZALPHA*ZALPHA)
         EBYM = 1.0D00/SQRT(1.0D00 + (ZALPHA/(GAMMA + NREL))**2)
         E(IORB) = (1.0D00 - EBYM)*CSQ
      ELSE
         E(IORB) = CSQ
      ENDIF
!
      IF (ZALPHA < FKABS) THEN
         GAMMA = SQRT(FKAP2 - ZALPHA*ZALPHA)
         EBYM = 1.0D00/SQRT(1.0D00 + (ZALPHA/(GAMMA + NREL - 0.5D00))**2)
         EMAX = (1.0D00 - EBYM)*CSQ
      ELSE
         EMAX = CSQ + CSQ
      ENDIF
!
      DELE = 0.0D00
!
!   Initialize
!
      FAIL = .FALSE.
      KOUNT = -1
!
!   Iteration loop begins here
!
    1 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT > MXK) THEN
         FAIL = .TRUE.
         RETURN
      ENDIF
!
!   Generate estimate of eigenvalue for this iteration
!
      EEST = E(IORB) + DELE
      IF (EEST>EMIN .AND. EEST<EMAX) THEN
         E(IORB) = EEST
      ELSE
         E(IORB) = 0.5D00*(EMIN + EMAX)
      ENDIF
!
!   Set up arrays TF and TG; find join point
!
      CALL SETPOT (IORB, JP)
!
!   Initialize outward integration
!
      ITYPE = 1
      CALL START (IORB, ITYPE, PZ(IORB), PF(:,IORB), Q0, QF(:,IORB))
!
!   Continue outward integration to classical turning point
!
      NSTRT = 6
      CALL SBSTEP (IORB, NSTRT, JP, PF(:,IORB), QF(:,IORB))
      PFJPO = PF(JP,IORB)
      QFJPO = QF(JP,IORB)
!
!   Initialize inward integration
!
      CALL TAIL (IORB, PF(:,IORB), QF(:,IORB), JP, MTP)
!
!   Continue inward integration to classical turning point
!
      NSTRT = MTP
      CALL SBSTEP (IORB, NSTRT, JP, PF(:,IORB), QF(:,IORB))
      PFJPI = PF(JP,IORB)
      QFJPI = QF(JP,IORB)
!
!   Make large component continuous, determine mismatch in small
!   component
!
      RATIO = PFJPO/PFJPI
      PF(JP:N,IORB) = PF(JP:N,IORB)*RATIO
      QF(JP:N,IORB) = QF(JP:N,IORB)*RATIO
!
!   Count nodes
!
      CALL COUNT (PF(:,IORB), MTP, NPC, SGN)
!
!   Correct the energy estimate if the number of nodes is wrong
!
      IF (NPC > NNP) THEN
         EMIN = E(IORB)
         DELE = 0.5D00*(EMAX - EMIN)
         GO TO 1
      ELSE IF (NPC < NNP) THEN
         EMAX = E(IORB)
         DELE = -0.5D00*(EMAX - EMIN)
         GO TO 1
      ENDIF
!
!   Correct number of nodes
!
!   Compute 'norm' of solution
!
      TA(1) = 0.0D00
      TA(2:N) = (PF(2:N,IORB)**2+QF(2:N,IORB)**2)*RP(2:N)
      MTP = N
      CALL QUAD (DNORM)
!
!   Determine correction to eigenvalue from magic formula
!   correct slope at origin
!
      QFJPI = QFJPI*RATIO
      DMSMCH = QFJPI - QFJPO
      IF (ABS(DMSMCH/QFJPO) > EPSLON) THEN
         DELE = C*PF(JP,IORB)*DMSMCH/DNORM
         IF (DELE < 0.0D00) THEN
            EMAX = EMAX*(1.0D00 - 0.2E00*ABS(DELE/E(IORB)))
         ELSE
            EMIN = EMIN*(1.0D00 + 0.2E00*ABS(DELE/E(IORB)))
         ENDIF
         PZ(IORB) = PZ(IORB)/SQRT(DNORM)
         GO TO 1
      ENDIF
!
!   Normalize
!
      DNFAC = 1.0D00/SQRT(DNORM)
      PZ(IORB) = PZ(IORB)*DNFAC
      PF(:N,IORB) = PF(:N,IORB)*DNFAC
      QF(:N,IORB) = QF(:N,IORB)*DNFAC
!
!   Find maximum tabulation point
!
      I = N + 1
    5 CONTINUE
      I = I - 1
      IF (ABS(PF(I,IORB)) < EPSLON) THEN
         PF(I,IORB) = 0.0D00
         QF(I,IORB) = 0.0D00
         GO TO 5
      ELSE
         MF(IORB) = I
      ENDIF
!
      RETURN
      END SUBROUTINE SOLVH
