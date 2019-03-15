!***********************************************************************
!                                                                      *
      SUBROUTINE BESSEL(IA, IB, IK, IW, K)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!                                                                      *
!   This routine evaluates the functions                               *
!                                                                      *
!                  (2K+1)!!                                            *
!      BESSJ  =    --------  J  (w   r) - 1  =  PHI  (w   r) - 1       *
!                         K   K   ab               K   ab              *
!                  (w   r)                                             *
!                    ab                                                *
!   and                                                                *
!                         K+1                                          *
!                  (w   r)                                             *
!                    ab                                                *
!      BESSN  =  - --------  N  (w   r) - 1  =  PSI  (w   r) - 1       *
!                  (2K-1)!!   K   ab               K   ab              *
!                                                                      *
!                                                                      *
!   where J and N ARE spherical Bessel functions, and                  *
!                                                                      *
!                w   = ABS ( E(IA) - E(IB) )/c                         *
!                 ab                                                   *
!                                                                      *
!   where  E(I)  is the eigenvalue for orbital  I . The writeup (B J   *
!   McKenzie, I P Grant, and P H Norrington, Computer Phys Commun 21   *
!   (1980) 233-246) is incorrect in its description of the output of   *
!   this routine.                                                      *
!                                                                      *
!   The  routine uses equations given in M Abramowitz and I A STegun   *
!   to evaluate the functions. Devices are used to reduce the number   *
!   of actual evaluations of these functions.                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M O D U L E S
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE bess_C
      USE debug_C
      USE def_C
      USE grid_C
      USE orb_C
      USE stor_C
      USE wfaC_c
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: IA
      INTEGER, INTENT(IN) :: IB
      INTEGER, INTENT(IN) :: IK
      INTEGER, INTENT(IN) :: IW
      INTEGER, INTENT(IN) :: K
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICODE, IWKP, IKKP, L, I, NN, J, JCHAN, IREM, ISWAP
      REAL(DOUBLE) :: EPSI, W, WA, XBESS1, XBESS2, S1, S2, DFNM, DFN, SSN, SCN&
         , SN, CN, OBWA, B, SKEEP
!-----------------------------------------------
!
      EPSI = DSQRT(0.1D00*ACCY)
!
!   Form unique label symmetric in IA, IB
!
      ICODE = MAX(IA,IB) + KEY*(MIN(IA,IB) + KEY*K)
!
!   Function in position; return
!
      IF (ICODE == KEEP(IK,IW)) RETURN
!
!   Function not in position; is it available in BESS arrays?
!
      W = WFACT*DABS(E(IA)-E(IB))/C
      WIJ(IW) = W
!
      DO IWKP = 1, 2
         DO IKKP = 1, 2
            IF (KEEP(IKKP,IWKP) /= ICODE) CYCLE
!
!   Function found move into position
!
            KEEP(IK,IW) = ICODE
            IF (LDBPR(7)) WRITE (99, 302) NP(IA), NH(IA), NP(IB), NH(IB), K, &
               IKKP, IWKP, IK, IW
            BESSJ(IK,IW,:N) = BESSJ(IKKP,IWKP,:N)
            BESSN(IK,IW,:N) = BESSN(IKKP,IWKP,:N)
            RETURN
         END DO
      END DO
!
!   Function not found; evaluate it
!
      IF (LDBPR(7)) WRITE (99, 303) NP(IA), NH(IA), NP(IB), NH(IB), K, IK, IW
!
      KEEP(IK,IW) = ICODE
!
      IF (W < EPSI**2) THEN
!
!   Negligible w
!
         BESSJ(IK,IW,:N) = 0.0D00
         BESSN(IK,IW,:N) = 0.0D00
         RETURN
!
      ENDIF
!
      NN = K
!
      BESSJ(IK,IW,1) = 0.0D00
      BESSN(IK,IW,1) = 0.0D00
!
!   Use a four-term power series for low w*r
!
      L5: DO J = 2, N
         WA = -0.5D00*(R(J)*W)**2
         XBESS1 = 1.0D00
         XBESS2 = 1.0D00
         S1 = 0.0D00
         S2 = 0.0D00
         DO I = 1, 4
            XBESS1 = XBESS1*WA/DBLE(I*(2*(NN + I) + 1))
            XBESS2 = XBESS2*WA/DBLE(I*(2*(I - NN) - 1))
            S1 = S1 + XBESS1
            S2 = S2 + XBESS2
            IF (DABS(XBESS1)>=DABS(S1)*EPSI .OR. DABS(XBESS2)>=DABS(S2)*EPSI) &
               CYCLE
            BESSJ(IK,IW,J) = S1
            BESSN(IK,IW,J) = S2
            CYCLE  L5
         END DO
         JCHAN = J
         GO TO 6
      END DO L5
!
!   If here then calculated whole array using four-term power
!   series.  Hence return
!
      RETURN
!
!   Use sin/cos expansion when power series requires more than
!   four terms terms to converge
!
    6 CONTINUE
      IF (NN == 0) THEN
         DFNM = 1.0D00
         DFN = 1.0D00
      ELSE
         DFNM = 1.0D00
         DO I = 3, 2*NN - 1, 2
            DFNM = DFNM*DBLE(I)
         END DO
         DFN = DFNM*DBLE(2*NN + 1)
      ENDIF
      DFNM = 1.0D00/DFNM
!
      IREM = MOD(NN,4)
!
      SELECT CASE (IREM)
      CASE (1)
!
!   NN = 1, 5, 9, ...
!
         SSN = -1.0D00
         SCN = 1.0D00
         ISWAP = 1
!
      CASE (2)
!
!   N = 2, 6, 10, ....
!
         SSN = -1.0D00
         SCN = -1.0D00
         ISWAP = 0
!
      CASE (3)
!
!   NN = 3, 7, 11,...
!
         SSN = 1.0D00
         SCN = -1.0D00
         ISWAP = 1
!
      CASE DEFAULT
!
!   NN = 0, 4, 8,...
!
         SSN = 1.0D00
         SCN = 1.0D00
         ISWAP = 0
!
      END SELECT
!
      DO J = JCHAN, N
         WA = W*R(J)
         IF (ISWAP == 0) THEN
            SN = SSN*DSIN(WA)
            CN = SCN*DCOS(WA)
         ELSE
            SN = SSN*DCOS(WA)
            CN = SCN*DSIN(WA)
         ENDIF
         OBWA = 1.0D00/WA
         B = OBWA
         S1 = B*SN
         S2 = B*CN
         DO I = 1, NN
            SKEEP = SN
            SN = CN
            CN = -SKEEP
            B = B*OBWA*DBLE((NN + I)*(NN - I + 1))/DBLE(2*I)
            S1 = S1 + B*SN
            S2 = S2 + B*CN
         END DO
         S1 = S1*DFN/WA**NN - 1.0D00
         S2 = S2*WA**(NN + 1)*DFNM - 1.0D00
         BESSJ(IK,IW,J) = S1
         BESSN(IK,IW,J) = S2
      END DO
      RETURN
!
  303 FORMAT(93X,I2,A2,2X,I2,A2,2X,I2,2X,'New',6X,'(',I2,',',I2,')')
  302 FORMAT(93X,I2,A2,2X,I2,A2,2X,I2,2X,'(',I2,',',I2,')',2X,'(',I2,',',I2,')'&
         )
      RETURN
!
      END SUBROUTINE BESSEL
