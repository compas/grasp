!***********************************************************************
!                                                                      *
      SUBROUTINE SOLVE(J, FAIL, INV, JP, NNP)
!                                                                      *
!   This subroutine  performs  step  2 in Algorithm 5.2 and 5.3 of C   *
!   Froese Fischer, Comput Phys Rep 3 (1986) 295. Some minor changes   *
!   have been made.                                                    *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      J     : (Input) The serial number of the orbital                *
!      JP    : (Output) The join point                                 *
!      FAIL  : (Output) If .true., the iterations did not yield an     *
!              acceptable solution (methods 1 and 2)                   *
!                                                                      *
!   Call(s) to: [RSCF92]: COUNT, ESTIM, IN, NEWE, OUT, cofpot,         *
!                         SETXUV, SETXV, SETXZ, START                  *
!               [LIB92]: DCBSRW, QUAD, SETPOT.                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
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
      USE debug_C
      USE def_C, ONLY: C, NSOLV
      USE grid_C, ONLY: rp
      USE int_C, ONLY: p, q, p0, q0, mtp0
      USE invt_C
      USE node_c
      USE orb_C
      USE scf_C
      USE tatb_C
      USE wave_C
      USE POTE_C
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE estim_I
      USE eigen_I
      USE dcbsrw_I
      USE setpot_I
      USE setxz_I
      USE start_I
      USE out_I
      USE in_I
      USE setxuv_I
      USE quad_I
      USE setxv_I
      USE prwf_I
      USE count_I
      USE newe_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      INTEGER  :: INV
      INTEGER  :: JP
      INTEGER  :: NNP
      LOGICAL  :: FAIL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KOUNT, NAKJ, NLOOPS, ICASE, MTPH, MTPI, I, MTPC, MTPV, MTPVC, &
         LOC1, LOC2, LOC, MX, NPRIME
      REAL(DOUBLE), DIMENSION(NNNP) :: PH, QH, PV, QV
      REAL(DOUBLE) :: TWOCSQ, ELAST, TENEJ, DELEPS, SGN, QJPOH, QJPIH, QJPOI, &
         QJPII, DNORM, ALFA, P0H, P0V, Q0V, QJPOV, QJPIV, CRNORM, DVNORM, AA, &
         BB, CC, DISCR, QQ, ROOT1, ROOT2, PMX, APM, PATI, ABPI, RATIO, TEST1, &
         TEST2, DELE
      LOGICAL :: CHECK
!-----------------------------------------------
!
!   Initialization
!
      CHECK = .NOT.NOINVT(J)
      FAIL = .FALSE.
      KOUNT = 0
      NAKJ = NAK(J)
!     TWOCSQ = 2.0D 00*137.036**2
      TWOCSQ = 2.0D00*C*C
!
!   Debug header
!
      IF (LDBPR(22)) WRITE (99, 300) NP(J), NH(J)
!
      NLOOPS = MAX(NSOLV,3*NP(J))
!
      CALL ESTIM (J)
      ELAST = E(J)
      E(J) = EIGEN(J)
!
!   Checks on lower bounds
!
      IF (E(J) < EPSMIN) THEN
         IF (E(J) > 0.0D00) THEN
            EPSMIN = E(J)
         ELSE
            WRITE (*, 301) NP(J), NH(J), E(J)
            E(J) = EPSMIN
            IF (ABS(EMIN - ELAST)<=1.0D-06 .AND. METHOD(J)<=2) THEN
               CALL DCBSRW (NP(J), NAKJ, ZINF, E(J), P0, P, Q, MTP0)
               WRITE (*, 302) NP(J), NH(J), ZINF
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
!   Check on upper bound
!
      IF (METHOD(J) <= 2) THEN
         IF (E(J) > EPSMAX) THEN
            TENEJ = 10.0D00*E(J)
            EPSMAX = MIN(TENEJ,TWOCSQ)
            EMAX = EPSMAX
            IF (E(J) > TWOCSQ) THEN
               WRITE (*, 303) NP(J), NH(J), E(J)
               E(J) = TWOCSQ
            ENDIF
         ENDIF
      ENDIF
!
!   Iteration loop begins here
!
    1 CONTINUE
      KOUNT = KOUNT + 1
!
      P0 = PZ(J)
!
      IF (KOUNT > 1) THEN
!
!   Check that bounds are ordered correctly
!
         IF (EPSMAX <= EPSMIN) WRITE (*, 304) EPSMIN, EPSMAX, NP(J), NH(J), E(J&
            )
!
         IF (KOUNT>NLOOPS .OR. EPSMAX-EPSMIN<1.0D00/DBLE(NP(J))**3) THEN
            WRITE (*, 305) METHOD(J), NP(J), NH(J)
            WRITE (*, 306) KOUNT - 1, NLOOPS, P0, E(J), DELEPS, EPSMIN, EPSMAX&
               , JP, MTP, NNP, NNODEP(J), SGN
            FAIL = .TRUE.
            RETURN
         ENDIF
      ENDIF
!
!   Set up arrays TF and TG; find join point
!
      CALL SETPOT (J, JP)
!
!   Set right-hand side to zero to form homogeneous equations;
!   integrate homogeneous equations outwards and inwards; store
!   small component at join point each time
!
      CALL SETXZ (J)
      ICASE = 1
      CALL START (J, ICASE, P0, PH, Q0, QH)
      CALL OUT (J, JP, PH, QH)
      QJPOH = QH(JP)
      CALL IN (J, JP, PH, QH, MTPH)
      QJPIH = QH(JP)
!
!   Set up right-hand side for inhomogeneous equations; integrate
!   inhomogeneous equations outwards and inwards; store small
!   component at join point each time
!
      CALL SETXUV (J)
      ICASE = 2
      CALL START (J, ICASE, P0, P, Q0, Q)
      CALL OUT (J, JP, P, Q)
      QJPOI = Q(JP)
      CALL IN (J, JP, P, Q, MTPI)
      QJPII = Q(JP)
!
!   Determine energy adjustment for methods 1 and 2
!
      IF (METHOD(J) <= 2) THEN
         TA(1) = 0.0D00
         TA(2:MTPI) = (P(2:MTPI)**2+Q(2:MTPI)**2)*RP(2:MTPI)
         MTP = MTPI
         CALL QUAD (DNORM)
!        DELEPS = 137.036*P(JP)*(QJPII-QJPOI)/DNORM
         DELEPS = C*P(JP)*(QJPII - QJPOI)/DNORM
      ENDIF
!
!   Generate the continuous solution
!
      MTPC = MAX(MTPH,MTPI)
      ALFA = -(QJPII - QJPOI)/(QJPIH - QJPOH)
      P0H = P0
      P0 = P0*(1.0D00 + ALFA)
      P(:MTPC) = P(:MTPC) + ALFA*PH(:MTPC)
      Q(:MTPC) = Q(:MTPC) + ALFA*QH(:MTPC)
!
      IF (METHOD(J)==2 .OR. METHOD(J)==4) THEN
!
!   Set up right-hand side for variational equations; integrate
!   variational equations outwards and inwards; store small
!   component at join point each time
!
         P0V = 0.0D00
         CALL SETXV (J)
         ICASE = 3
         CALL START (J, ICASE, P0V, PV, Q0V, QV)
         CALL OUT (J, JP, PV, QV)
         QJPOV = QV(JP)
         CALL IN (J, JP, PV, QV, MTPV)
         QJPIV = QV(JP)
!
!   Generate continuous solutions
!
         MTPVC = MAX(MTPC,MTPV)
         ALFA = -(QJPIV - QJPOV)/(QJPIH - QJPOH)
         PV(:MTPVC) = PV(:MTPVC) + ALFA*PH(:MTPVC)
         QV(:MTPVC) = QV(:MTPVC) + ALFA*QH(:MTPVC)
!
         TA(1) = 0.0D00
         TA(2:MTPC) = RP(2:MTPC)*(P(2:MTPC)**2+Q(2:MTPC)**2)
         MTP = MTPC
         CALL QUAD (DNORM)
!
         MTP = MIN(MTPC,MTPVC)
         TA(1) = 0.0D00
         TA(2:MTP) = RP(2:MTP)*(P(2:MTP)*PV(2:MTP)+Q(2:MTP)*QV(2:MTP))
         CALL QUAD (CRNORM)
!
         TA(1) = 0.0D00
         TA(2:MTPVC) = RP(2:MTPVC)*(PV(2:MTPVC)**2+QV(2:MTPVC)**2)
         MTP = MTPVC
         CALL QUAD (DVNORM)
!
!   Determine deleps required to normalize new solution to
!   first order: modified form of solution to a quadratic
!   equation (see Press et al.)
!
         AA = DVNORM
         BB = CRNORM + CRNORM
         CC = DNORM - 1.0D00
         DISCR = BB*BB - 4.0D00*AA*CC
         IF (DISCR > 0.0D00) THEN
            QQ = -0.5D00*(BB + SIGN(1.0D00,BB)*SQRT(DISCR))
            ROOT1 = CC/QQ
            ROOT2 = QQ/AA
            PMX = 0.0D00
            APM = 0.0D00
            DO I = 2, JP
               PATI = P(I)
               IF (PATI > PMX) THEN
                  PMX = PATI
                  LOC1 = I
               ENDIF
               ABPI = ABS(PATI)
               IF (ABPI <= APM) CYCLE
               APM = ABPI
               LOC2 = I
            END DO
            IF (PMX /= 0.0D00) THEN
               RATIO = APM/ABS(PMX)
               IF (RATIO < 10.0D00) THEN
                  LOC = LOC1
               ELSE
                  LOC = LOC2
               ENDIF
            ELSE
               LOC = LOC2
            ENDIF
            TEST1 = P(LOC) + ROOT1*PV(LOC)
            TEST2 = P(LOC) + ROOT2*PV(LOC)
            IF (TEST1>0.0D00 .AND. TEST2<0.0D00) THEN
               DELE = ROOT1
            ELSE IF (TEST1<0.0D00 .AND. TEST2>0.0D00) THEN
               DELE = ROOT2
            ELSE IF (TEST1>0.0D00 .AND. TEST2>0.0D00) THEN
               IF (TEST1 < TEST2) THEN
                  DELE = ROOT1
               ELSE
                  DELE = ROOT2
               ENDIF
            ELSE IF (TEST1<0.0D00 .AND. TEST2<0.0D00) THEN
               IF (TEST1 > TEST2) THEN
                  DELE = ROOT1
               ELSE
                  DELE = ROOT2
               ENDIF
            ENDIF
         ELSE
            DELE = -BB/(AA + AA)
         ENDIF
!
!   Generate new solution
!
         MTP0 = MAX(MTPC,MTPVC)
         P0 = P0 + DELE*ALFA*P0H
         P(2:MTP0) = P(2:MTP0) + DELE*PV(2:MTP0)
         Q(2:MTP0) = Q(2:MTP0) + DELE*QV(2:MTP0)
      ELSE
         MTP0 = MTPC
      ENDIF
!
!   Debug printout
!
      IF (LDBPR(23)) CALL PRWF (J)
!
!   Count nodes in large component function; determine sign at first
!   oscillation, effective quantum number; note that node counting
!   is never enforced on the small component
!
      CALL COUNT (P, MTP0, NNP, SGN)
!
!   DEBUG PRINTOUT
!
      IF (LDBPR(22)) WRITE (99, 306) KOUNT, NLOOPS, P0, E(J), DELEPS, EPSMIN, &
         EPSMAX, JP, MTP, NNP, NNODEP(J), SGN
!
!   Proceed according to method
!
      IF (METHOD(J) > 2) THEN
         IF (CHECK .AND. SGN<0.0D00) THEN
            INV = 1
            P0 = -P0
            P(2:MTP0) = -P(2:MTP0)
            Q(2:MTP0) = -Q(2:MTP0)
         ENDIF
      ELSE
         MX = NNP - NNODEP(J)
         NPRIME = NNP + NKL(J) + 1
         CALL NEWE (J, SGN, NPRIME, MX, DELEPS, FAIL, INV)
         IF (FAIL) GO TO 1
      ENDIF
!
!   Solution found
!
      RETURN
!
  300 FORMAT(/,/,' Debug printout active; orbital: ',1I2,1A2)
  301 FORMAT(' E(',1I2,1A2,') = ',1P,D11.4,'; adjusted to EPSMIN')
  302 FORMAT(' Returned hydrogenic function for ',1I2,1A2,' with',&
         ' effective charge ',F7.3)
  303 FORMAT(' E(',1I2,1A2,') = ',1P,D11.4,'; adjusted to TWOCSQ')
  304 FORMAT(' Warning: difficulty with node-counting procedure'/,&
         ' lower bound on energy (',1P,D11.4,') exceeds upper',' bound (',1D&
         11.4,'; E(',1I2,1A2,') = ',1D11.4)
  305 FORMAT(' Method ',1I1,' unable to solve for ',1I2,1A2,' orbital')
  306 FORMAT(' Iteration number: ',1I2,', limit: ',1I2,/,&
         ' Present estimate of P0; ',1D21.14,/,' Present estimate of E(J): ',1D&
         21.14,', DELEPS: ',1D21.14,/,' Lower bound on energy: ',1D21.14,&
         ', upper bound: ',1D21.14,/,' Join point: ',1I4,&
         ', Maximum tabulation point:',1I4,/,' Number of nodes counted: ',1I2,&
         ', Correct number: ',1I2,/,' Sign of P at first oscillation: ',F3.0)
      RETURN
!
      END SUBROUTINE SOLVE
