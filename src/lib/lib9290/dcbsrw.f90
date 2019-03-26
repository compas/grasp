!***********************************************************************
!                                                                      *
      SUBROUTINE DCBSRW(N, KAPPA, Z, E, RG0, RG, RF, MTP)
!                                                                      *
!   This subroutine computes the  Dirac-Coulomb  bound-state orbital   *
!   radial wavefunction.   Equations (13.5) and (13.5') of  Akhiezer   *
!   and Berestetskii modified to ensure positive slope at the origin   *
!   for RG are used.                                                   *
!                                                                      *
!   The arguments are as follows:                                      *
!                                                                      *
!      N    : (Input) The (usual) principal quantum number             *
!      KAPPA: (Input) The relativistic angular quantum number          *
!      Z    : (Input) The effective nuclear charge                     *
!      E    : (Output) The Dirac-Coulomb Eigenenergy                   *
!      RG0  : (Output) Coefficient of the leading term in the          *
!                      series expansion of the large component         *
!                      near the origin                                 *
!      RG   : (Output) r times the large component wavefunction of     *
!                      Akhiezer and Berestetskii                       *
!      RF   : (Output) r times the small component wavefunction of     *
!                      Akhiezer and Berestetskii                       *
!      MTP  : (Output) Maximum tabulation point                        *
!                                                                      *
!   Call(s) to: [LIB92]: CGAMMA.                                       *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def, ONLY: NNNP
      USE DEF_C ,          ONLY: C, ACCY
      USE GRID_C, ONLY: R, NTP=>N
      USE TATB_C, ONLY: TA, TB
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cgamma_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: KAPPA
      INTEGER, INTENT(OUT) :: MTP
      REAL(DOUBLE), INTENT(IN) :: Z
      REAL(DOUBLE), INTENT(OUT) :: E
      REAL(DOUBLE), INTENT(OUT) :: RG0
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(INOUT) :: RG
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(OUT) :: RF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, NR, NRFAC, I, IORDR1, IORDR2
      REAL(DOUBLE) :: ALFA, FN, FKAPPA, FK, FNR, ZALFA, GAMMA, TWOGP1, BIGN, &
         EPS, ARGR, ARGI, RGAMM1, DUMMY, RGAMM2, FAC, FG, FF, FACN, A, AN1, AN2&
         , B, BN, FDEN, BIGNMK, RHO, RHON, F1, F2, OVLFAC, CUTOFF
!-----------------------------------------------
!
!
!
!   Ensure that the principal quantum number is physical
!
      IF (N <= 0) THEN
         WRITE (*, 300)
         WRITE (*, 301) N
         STOP
      ENDIF
!
!   Ensure that the angular quantum number is physical
!
      IF (KAPPA == 0) THEN
         WRITE (*, 300)
         WRITE (*, 302)
         STOP
      ELSE IF (KAPPA == N) THEN
         WRITE (*, 300)
         WRITE (*, 303) KAPPA, N
         STOP
      ELSE IF (ABS(KAPPA) > N) THEN
         WRITE (*, 300)
         WRITE (*, 303) KAPPA, N
         STOP
      ENDIF
!
!   Ensure that the charge is physical
!
      IF (Z <= 0.0D00) THEN
         WRITE (*, 300)
         WRITE (*, 304) Z
         STOP
      ELSE IF (Z > C) THEN
         WRITE (*, 300)
         WRITE (*, 305) Z, C
         STOP
      ENDIF
!
!   Atomic units
!
      ALFA = 1.0D00/C
!
!   Now determine all the parameters
!
      FN = DBLE(N)
      FKAPPA = DBLE(KAPPA)
      K = ABS(KAPPA)
      FK = DBLE(K)
      NR = N - K
      FNR = DBLE(NR)
      ZALFA = Z*ALFA
      GAMMA = SQRT(FK*FK - ZALFA*ZALFA)
      TWOGP1 = GAMMA + GAMMA + 1.0D00
      BIGN = SQRT(FN*FN - 2.0D00*FNR*(FK - GAMMA))
      EPS = 1.0D00/SQRT(1.0D00 + (ZALFA/(GAMMA + FNR))**2)
!
!   EPS is the total energy divided by C*C; this must be converted
!   to the units and reference energy of GRASP
!
      E = (1.0D00 - EPS)*C*C
!
!   Now the normalization constants
!
      NRFAC = 1
      DO I = 1, NR
         NRFAC = NRFAC*I
      END DO
!
      ARGR = TWOGP1 + FNR
      ARGI = 0.0D00
      CALL CGAMMA (ARGR, ARGI, RGAMM1, DUMMY)
      ARGR = TWOGP1
      CALL CGAMMA (ARGR, ARGI, RGAMM2, DUMMY)
!
      FAC = -SQRT(RGAMM1)/(RGAMM2*SQRT(DBLE(NRFAC)))*SQRT(Z/(2.0D00*BIGN*BIGN*(&
         BIGN - FKAPPA)))
!
!   Ensure that the slope of the large-component function is
!   positive at the origin
!
      IF (KAPPA > 0) FAC = -FAC
!
      FG = FAC*SQRT(1.0D00 + EPS)
      FF = FAC*SQRT(1.0D00 - EPS)
!
!   Now set up the coefficients of the confluent hypergeometric
!   functions  F (-NR+1,2*GAMMA+1;RHO)  and  F (-NR,2*GAMMA+1;RHO)
!   in the workspace arrays  TA  and  TB , respectively
!
      IF (NR == 0) THEN
         IORDR1 = 0
         IORDR2 = 0
      ELSE
         IORDR1 = NR - 1
         IORDR2 = NR
      ENDIF
!
      FAC = 1.0D00
      FACN = 1.0D00
      A = -FNR
      AN1 = A + 1.0D00
      AN2 = A
      B = TWOGP1
      BN = B
!
      K = 0
    2 CONTINUE
      K = K + 1
      FDEN = 1.0D00/(FACN*BN)
      IF (K <= IORDR1) TA(K) = AN1*FDEN
      IF (K <= IORDR2) THEN
         TB(K) = AN2*FDEN
         A = A + 1.0D00
         AN1 = AN1*(A + 1.0D00)
         AN2 = AN2*A
         B = B + 1.0D00
         BN = BN*B
         FAC = FAC + 1.0D00
         FACN = FACN*FAC
         GO TO 2
      ENDIF
!
!   Now tabulate the function over the entire grid
!
      RG(1) = 0.0D00
      RF(1) = 0.0D00
      FAC = (Z + Z)/BIGN
      BIGNMK = BIGN - FKAPPA
      DO I = 2, NTP
         RHO = FAC*R(I)
         RHON = RHO
         K = 0
         F1 = 1.0D00
         F2 = 1.0D00
    3    CONTINUE
         K = K + 1
         IF (K <= IORDR1) F1 = F1 + TA(K)*RHON
         IF (K <= IORDR2) THEN
            F2 = F2 + TB(K)*RHON
            RHON = RHON*RHO
            GO TO 3
         ENDIF
         F1 = FNR*F1
         F2 = BIGNMK*F2
         OVLFAC = EXP((-0.5D00*RHO))*RHO**GAMMA
         RG(I) = FG*OVLFAC*(F1 - F2)
         RF(I) = FF*OVLFAC*(F1 + F2)
      END DO
!
!   Determine the effective maximum tabulation point based on the
!   cutoff; define the cutoff conservatively
!
      CUTOFF = ACCY*0.1D00
!
      MTP = NTP + 1
    5 CONTINUE
      MTP = MTP - 1
!      IF (ABS(RG(MTP)) < CUTOFF) THEN
      IF ((ABS(RG(MTP)) < CUTOFF).OR.(ABS(R(MTP)) > 1.D+50)) THEN ! JE: APPLY BOX (R<1.D+50)
         RG(MTP) = 0.0D00
         RF(MTP) = 0.0D00
         GO TO 5
      ENDIF
!
      IF (MTP == NTP) WRITE (*, 306) NTP, RG(NTP), CUTOFF
!
!   Compute the coefficient of R**GAMMA at the origin
!
      RG0 = FG*FAC**GAMMA*(FNR - BIGNMK)
!
      RETURN
!
  300 FORMAT('DCBSRW:')
  301 FORMAT(' Principal quantum number is ',1I3)
  302 FORMAT(' Angular quantum number is 0')
  303 FORMAT(' Angular quantum number (',1I3,') is out of range for',&
         ' principal quantum number (',1I3,')')
  304 FORMAT(' Nuclear charge (',3P,1D16.7,') is too small')
  305 FORMAT(' Nuclear charge (',3P,1D16.7,') exceeds C (',1D16.7,')')
  306 FORMAT(/,/,/,' ***** Warning in SUBROUTINE DCBSRW *****'/,/,&
         ' Radial grid of insufficient extent:'/,' P(',1I4,') = ',1P,1D10.3,&
         ', Exceeds cutoff (',1D10.3,')')
      RETURN
!
      END SUBROUTINE DCBSRW
