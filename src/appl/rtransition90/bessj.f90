!***********************************************************************
!                                                                      *
      SUBROUTINE BESSJ(W)
!                                                                      *
!   This routine evaluates Bessel fuctions J K ( W*R/C ) at the grid   *
!   points for K=L-1,L,L+1 and stores  them in the  arrays BJ(..,1),   *
!   BJ(..,2),BJ(..,3) respectively. It uses a power series expansion   *
!   for  small  r and switches to sin/cos expansion when more than 4   *
!   terms in power series are required.                                *
!                                         Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE bess_C, ONLY: BJ, TC, TD
      USE debug_C, ONLY: LDBPR
      USE grid_C
      USE osc_C, ONLY: LK, KK, L=>LK
      USE DEF_C, ONLY: C, CVAC, PI
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN) :: W
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I1, IW, NN, IEND, IPROD, I, JCHAN, J, ISWAP, MODNN4, &
                 JJ, JJJ, II
      REAL(DOUBLE) :: EPSI, DFN, WR, WA, XBESS1, S1, SSN, SCN, SN, CN, B, SKEEP
!-----------------------------------------------
!
!
      EPSI = 1.0D-05
!
      IW = 1
      NN = L - 1
      IF (KK /= 0) THEN
         IW = 2
         NN = L
      ENDIF
    1 CONTINUE
      IEND = 2*NN + 1
      IPROD = 1
      DO I = 1, IEND, 2
         IPROD = IPROD*I
      END DO
      DFN = IPROD
      JCHAN = N
      BJ(1,IW) = 1.0D00
      DO J = 2, N
         WR = W*R(J)
         WA = -WR*WR*0.5D00
         XBESS1 = 1.0D00
         S1 = 1.0D00
         DO I = 1, 4
            XBESS1 = XBESS1*WA/DBLE(I*(2*(NN + I) + 1))
            S1 = S1 + XBESS1
            IF (ABS(XBESS1) < ABS(S1)*EPSI) GO TO 4
         END DO
         JCHAN = J
         EXIT
    4    CONTINUE
         BJ(J,IW) = S1*WR**NN/DFN
      END DO
!
!   Use sin/cos expansion when power series takes longer
!   than 4 terms to converge
!
      IF (JCHAN < N) THEN
         ISWAP = 0
         MODNN4 = MOD(NN - 1,4) + 1
         SELECT CASE (MODNN4)
!
!   NN = 1, 5, 9, ...
!
         CASE DEFAULT
            SSN = -1.0D00
            SCN = 1.0D00
            ISWAP = 1
!
!   N = 2, 6, 10,....
!
         CASE (2)
            SSN = -1.0D00
            SCN = -1.0D00
!
!   NN = 3, 7, 11,...
!
         CASE (3)
            SSN = 1.0D00
            SCN = -1.0D00
            ISWAP = 1
!
!   NN = 0, 4, 8,...
!
         CASE (4)
            SSN = 1.0D00
            SCN = 1.0D00
         END SELECT
!
   13    CONTINUE
         DO J = JCHAN, N
            WA = W*R(J)
            IF (ISWAP <= 0) THEN
               SN = SSN*SIN(WA)
               CN = SCN*COS(WA)
            ELSE
               SN = SSN*COS(WA)
               CN = SCN*SIN(WA)
            ENDIF
            I = -1
            S1 = 0.0D00
            I = I + 1
            I1 = I
            DO I = I1, NN
               IF (I == 0) THEN
                  B = 1.0D00/WA
               ELSE
                  B = B*DBLE((NN + I)*(NN - I + 1))/DBLE(2*I)/WA
               ENDIF
               S1 = S1 + B*SN
               SKEEP = SN
               SN = CN
               CN = -SKEEP
            END DO
            BJ(J,IW) = S1
         END DO
      ENDIF
      IF (NN>=L + 1 .OR. KK==1) THEN
!
!   Print out Bessel functions if (debug) option set
!
         IF (LDBPR(16)) THEN
            DO JJ = 1, 3
               JJJ = L - 2 + JJ
               WRITE (99, 300) JJJ, (BJ(II,JJ),II=1,N)
            END DO
         ENDIF
!
!   All done
!
! zou
         DO JJ = 1, 3
            JJJ = L - 2 + JJ
            BJ(:N,JJ) = BJ(:N,JJ)*(C/CVAC)**JJJ
         END DO
! zou
         RETURN
      ELSE
         NN = NN + 1
         IW = IW + 1
         GO TO 1
      ENDIF
!
  300 FORMAT(/,' Bessel function of order ',I3,/(1P,7D18.10))
      RETURN
!
      END SUBROUTINE BESSJ
