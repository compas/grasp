

!***********************************************************************
!                                                                      *
      SUBROUTINE TAIL(IORB, P, Q, JP, MTP)
!                                                                      *
!   This subroutine begins the inward integration of the homogeneous   *
!   Dirac radial equation. With only minor modifications, the series   *
!   given by  J E Sienkiewicz  and W E Baylis, J Phys B: At Mol Phys   *
!   20 (1987) 5145-5156, p 5155, is used.                              *
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
      USE parameter_def, ONLY: NNNP
      USE DEF_C, ONLY: ACCY, C, EXPMAX
      USE GRID_C
      USE ORB_C
      USE POTE_C, ONLY: YP
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IORB
      INTEGER , INTENT(IN) :: JP
      INTEGER , INTENT(OUT) :: MTP
      REAL(DOUBLE) , INTENT(INOUT) :: P(NNNP)
      REAL(DOUBLE) , INTENT(OUT) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NM4, M, LOC
      REAL(DOUBLE) :: EPS, BIGE, BIGEBC, BGEBC2, BETA, FK, RJP, QEM, OFFSET, T&
         , ZLOC, GAM2, ZLBB, DNU, FKMZBB, SUMP, SUMQ, EM, CM, EMFACT, OVLTRM, &
         PTERM, QTERM, EXPTRM
!-----------------------------------------------
!
!
!   Initialize
!
      EPS = ACCY*0.1D00
      BIGE = -E(IORB)
      BIGEBC = BIGE/C
      BGEBC2 = BIGEBC/C
      BETA = SQRT((-BIGE*(2.0D00 + BGEBC2)))
      FK = DBLE(NAK(IORB))
!
!   Find MTP
!
      I = JP
      NM4 = N - 4
      RJP = R(JP)
      QEM = 0.25D00*EXPMAX
    1 CONTINUE
      I = I + 1
      IF (I <= NM4) THEN
         IF (BETA*(R(I)-RJP) > QEM) THEN
            MTP = I
         ELSE
            GO TO 1
         ENDIF
      ELSE
         WRITE (*, 300)
         MTP = NM4
      ENDIF
!
!   Compute offset for exponential function
!
      OFFSET = BETA*R(MTP)
!
!   Tabulate tail points
!
      DO I = MTP, N
!
         T = 2.0D00*BETA*R(I)
!
         ZLOC = YP(I)
         GAM2 = FK**2 - (ZLOC/C)**2
         ZLBB = ZLOC/BETA
         DNU = ZLBB*(1.0D00 + BGEBC2)
         FKMZBB = FK - ZLBB
!
         M = -1
         SUMP = 0.0D00
         SUMQ = 0.0D00
         M = M + 1
         EM = DBLE(M)
         IF (M == 0) THEN
            CM = 1.0D00
            EMFACT = 1.0D00
         ELSE
            CM = CM*(GAM2 - (DNU - EM)**2)
            EMFACT = EMFACT*EM
         ENDIF
         OVLTRM = CM*(T**(DNU - EM)/EMFACT)
         PTERM = OVLTRM*(FKMZBB + EM)*BETA
         QTERM = OVLTRM*(FKMZBB - EM)*BIGEBC
         SUMP = SUMP + PTERM
         SUMQ = SUMQ + QTERM
         DO WHILE(ABS(PTERM/SUMP)>=EPS .OR. ABS(QTERM/SUMQ)>=EPS)
            M = M + 1
            EM = DBLE(M)
            IF (M == 0) THEN
               CM = 1.0D00
               EMFACT = 1.0D00
            ELSE
               CM = CM*(GAM2 - (DNU - EM)**2)
               EMFACT = EMFACT*EM
            ENDIF
            OVLTRM = CM*(T**(DNU - EM)/EMFACT)
            PTERM = OVLTRM*(FKMZBB + EM)*BETA
            QTERM = OVLTRM*(FKMZBB - EM)*BIGEBC
            SUMP = SUMP + PTERM
            SUMQ = SUMQ + QTERM
         END DO
         EXPTRM = EXP((-0.5D00*T) + OFFSET)
         P(I) = SUMP*EXPTRM
         Q(I) = SUMQ*EXPTRM
         IF (P(I) /= 0.0D00) CYCLE
         LOC = I + 1
         GO TO 4
      END DO
      LOC = N + 1
!
    4 CONTINUE
      P(LOC:N) = 0.0D00
      Q(LOC:N) = 0.0D00
!
      RETURN
!
  300 FORMAT('TAIL: Grid may be of insufficient extent')
      RETURN
!
      END SUBROUTINE TAIL
