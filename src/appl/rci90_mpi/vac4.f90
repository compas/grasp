!***********************************************************************
!                                                                      *
      SUBROUTINE VAC4
!                                                                      *
!   This routine sets up the fourth-order vacuum polarization poten-   *
!   tial using equations (11) and (12) of L Wayne Fullerton and  G A   *
!   Rinker, Jr,  Phys  Rev  A 13 (1976) 1283-1287. The potential  is   *
!   accumulated in array  TC(I), I = 1, ..., N. It is transferred to   *
!   array TA in COMMON block TATB.                                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: FUNL.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE debug_C
      USE def_C
      USE grid_C, ONLY: r, n
      USE npar_C
      USE ncdist_C, ONLY: zdist
      USE tatb_C, ONLY: mtp, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE funl_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K, NB2, NROWS, II, II1, II2
      REAL(DOUBLE), DIMENSION(NNNP) :: TC
      REAL(DOUBLE) :: EPSI, TWOCV, FACTOR, RI, X, TCI, RK, XK, XI, XM, XP
!-----------------------------------------------
!
!   Overall initialization
!
      EPSI = PRECIS*PRECIS
      TWOCV = CVAC + CVAC
!
!   Potential for point nucleus: equation (12)
!
      FACTOR = -Z/(PI*CVAC)**2
!
      TC(1) = 0.0D00
!
      I = 1
    1 CONTINUE
      I = I + 1
      RI = R(I)
      X = TWOCV*RI
      TCI = (FACTOR/RI)*FUNL(X,1)
      IF (DABS(TCI) >= EPSI) THEN
         TC(I) = TCI
         IF (I < N) GO TO 1
      ELSE
         TC(I:N) = 0.0D00
      ENDIF
!
!   Potential for finite nucleus: equation (11)
!
      IF (NPARM == 2) THEN
!
         FACTOR = -1.0D00/(PI*CVAC**3)
!
         TC(1) = 0.0D00
!
         K = 1
    3    CONTINUE
         K = K + 1
!
         RK = R(K)
         XK = TWOCV*RK
         TA(1) = 0.0D00
!
         DO I = 2, MTP
            XI = TWOCV*R(I)
            XM = DABS(XK - XI)
            XP = XK + XI
            TA(I) = (FUNL(XM,0) - FUNL(XP,0))*ZDIST(I)
         END DO
!
         CALL QUAD (X)
!
         X = X*FACTOR/RK
!
!   Get out of the loop if the asymptotic region has been reached
!
         IF (DABS(X) >= EPSI) THEN
            IF (DABS((TC(K)-X)/X) > 1.0D-03) THEN
               TC(K) = X
               IF (K < N) GO TO 3
            ENDIF
         ENDIF
!
      ENDIF
!
      IF (LDBPR(8)) THEN
         WRITE (99, 300)
         NB2 = N/2
         IF (2*NB2 == N) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2 + 1
         ENDIF
         DO II = 1, NROWS
            II1 = II
            II2 = II1 + NROWS
            IF (II2 <= N) THEN
               WRITE (99, 301) R(II1), TB(II1), TC(II1), R(II2), TB(II2), TC(&
                  II2)
            ELSE IF (II1 <= N) THEN
               WRITE (99, 301) R(II1), TB(II1), TC(II1)
            ENDIF
         END DO
      ENDIF
!
!   Generate total vacuum-polarization potential
!
      TB(:N) = TC(:N) + TB(:N)
!
      RETURN
!
  300 FORMAT(/,/,/,' ++++++++++ VAC4 ++++++++++'/,/,2(&
         ' -------- r -------- ----- VV2 (r) -----',' ----- VV4 (r) -----'))
  301 FORMAT(1P,6(1X,1D19.12))
      RETURN
!
      END SUBROUTINE VAC4
