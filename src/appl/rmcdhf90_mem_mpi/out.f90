!***********************************************************************
!                                                                      *
      SUBROUTINE OUT(J, JP, P, Q)
!                                                                      *
!   This subroutine carries out the step-by-step outward integration   *
!   of a pair of inhomogeneous Dirac radial equations.                 *
!                                                                      *
!   arguments:                                                         *
!                                                                      *
!      J:   (Input) Orbital index of function to be computed           *
!      JP:  (Input) The join point; the outward integration stops      *
!           at this tabulation index                                   *
!      P,Q: (Input and output) on input, elements 1 to 3 of both       *
!           arrays must be tabulated; on output, the arrays are        *
!           tabulated up to point JP                                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP
      USE grid_C,           ONLY: h, rpor
      USE int_C,            ONLY: TF ,TG ,XU, XV
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
      INTEGER , INTENT(IN) :: JP
      REAL(DOUBLE) , INTENT(INOUT) :: P(NNNP)
      REAL(DOUBLE) , INTENT(INOUT) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: DKHB2, DKHB2F, CPI, CMI, PI, QI, TFI, TGI, CPIM1, CMIM1, &
         PIM1, QIM1, TFIM1, TGIM1, UCIM1, UDIM1, UEI
!-----------------------------------------------
!
!   One global initialization
!
      DKHB2 = 0.5D00*H*DBLE(NAK(J))
!
!   Tabulate P(r) and Q(r) by step-by-step integration
!
!   Initializations: set quantities for I = 3
!
      I = 3
      DKHB2F = DKHB2*RPOR(I)
      CPI = 1.0D00 + DKHB2F
      CMI = 1.0D00 - DKHB2F
      PI = P(I)
      QI = Q(I)
      TFI = TF(I)
      TGI = TG(I)
!
!   March out to from I = 4 to I = JP
!
!XHH Use doo-loop
      DO I = 4, JP
         CPIM1 = CPI
         CMIM1 = CMI
         DKHB2F = DKHB2*RPOR(I)
         CPI = 1.0D00 + DKHB2F
         CMI = 1.0D00 - DKHB2F
         PIM1 = PI
         QIM1 = QI
         TFIM1 = TFI
         TGIM1 = TGI
         TFI = TF(I)
         TGI = TG(I)
         UCIM1 = CMIM1*PIM1 - TFIM1*QIM1 + XU(I-1)
         UDIM1 = CPIM1*QIM1 - TGIM1*PIM1 + XV(I-1)
         UEI = CPI*CMI - TFI*TGI
         PI = (CMI*UCIM1 - TFI*UDIM1)/UEI
         QI = (CPI*UDIM1 - TGI*UCIM1)/UEI
         P(I) = PI
         Q(I) = QI
      END DO
!
      RETURN
      END SUBROUTINE OUT
