!***********************************************************************
!                                                                      *
      SUBROUTINE VAC2
!                                                                      *
!   This routine sets up the second-order vacuum polarization poten-   *
!   tial using  equations (1) and (4) of  L Wayne Fullerton and  G A   *
!   Rinker, Jr,  Phys Rev A  13 (1976) 1283-1287.  The  potential is   *
!   accumulated  in  array  TB(I), I = 1, ..., N  which is in COMMON   *
!   block /TATB/ .                                                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: FUNK.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE grid_C
      USE npar_C
      USE ncdist_C
      USE tatb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE funk_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K
      REAL(DOUBLE) :: EPSI, TWOCV, FACTOR, RI, X, TBI, RK, XK, XI, XM, XP
!-----------------------------------------------
!
!   Overall initialization
!
      EPSI = PRECIS*PRECIS
      TWOCV = CVAC + CVAC
!
!   Potential for a point nucleus: equation (1)
!   (this is also the asymptotoc form for a finite nucleus)
!
      FACTOR = -(2.0D00*Z)/(3.0D00*PI*CVAC)
!
      TB(1) = 0.0D00
!
      I = 1
    1 CONTINUE
      I = I + 1
!
      RI = R(I)
      X = TWOCV*RI
      TBI = (FACTOR/RI)*FUNK(X,1)
!
      IF (DABS(TBI) >= EPSI) THEN
         TB(I) = TBI
         IF (I < N) GO TO 1
      ELSE
         TB(I:N) = 0.0D00
      ENDIF
!
!   Potential for a finite nucleus: equation (4)
!
      IF (NPARM == 2) THEN
!
         FACTOR = -2.0D00/(3.0D00*CVAC**2)
!
!   Set up integrand
!
         TB(1) = 0.0D00
!
         K = 1
    3    CONTINUE
         K = K + 1
!
         RK = R(K)
         XK = TWOCV*RK
!
         TA(1) = 0.0D00
         DO I = 2, MTP
            XI = TWOCV*R(I)
            XM = DABS(XK - XI)
            XP = XK + XI
            TA(I) = (FUNK(XM,0) - FUNK(XP,0))*ZDIST(I)
         END DO
!
         CALL QUAD (X)
!
         X = X*FACTOR/RK
!
!   Get out of loop if the asymptotic value has been attained
!
         IF (DABS(X) >= EPSI) THEN
            IF (DABS((X - TB(K))/X) > 1.0D-05) THEN
               TB(K) = X
               IF (K < N) GO TO 3
            ENDIF
         ENDIF
!
      ENDIF
!
      RETURN
      END SUBROUTINE VAC2
