!***********************************************************************
      SUBROUTINE IMPROV(EOL, J, LSORT, DAMPMX) 
!************************************************************************     
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CORRE_C 
      USE damp_C
      USE def_C
      USE grid_C
      USE int_C
      USE mpi_s
      USE orb_C
      USE orthct_C
      USE scf_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cofpot_I 
      USE defcor_I 
      USE solve_I 
      USE orthsc_I 
      USE matrix_I 
      USE newco_I 
      USE setlag_I 
      USE quad_I 
      USE consis_I 
      USE dampck_I 
      USE dampor_I 
      USE orthy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
      REAL(DOUBLE) , INTENT(INOUT) :: DAMPMX 
      LOGICAL  :: EOL, LSORT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: P2 = 2.0D-01 
      REAL(DOUBLE), PARAMETER :: P005 = 5.0D-03 
      REAL(DOUBLE), PARAMETER :: P0001 = 1.0D-04 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IPR, NPTS, INV, JP, NNP, I, NWWW 
      REAL(DOUBLE) :: ED1, GAMAJ, ED2, WTAEV, DNORM, DNFAC, DEL1, DEL2, ODAMPJ 
      LOGICAL :: FAIL, FIRST 
!-----------------------------------------------
!
!   C Froese Fischer's IPR and ED1 parameter
!
      DATA IPR/ 0/  
      DATA ED1/ 0.D0/  
      DATA FIRST/ .FALSE./  
!
!
!-----------------------------------------------------------------------
      GAMAJ = GAMA(J) 
!
!   C Froese Fischer's parameters IPR, ED1, ED2 are set and
!   used in this routine and in DAMPCK
!
    1 CONTINUE 
      ED2 = E(J) 
      ED1 = PED(J)
!
!   Set up the exchange potential and arrays XU, XV as appropriate
!
!   Set coefficients for YPOT, XPOT, DACON
!   Compute direct potential, exchange potential
!   Add in Lagrange-multiplier contribution
!   Add in derivative-terms contribution
!
      NPTS = N 
      CALL COFPOT (EOL, J, NPTS) 
!
!   Calculate deferred corrections
!
      CALL DEFCOR (J) 
!
!   Solve the Dirac equation
!
      INV = 0 
      CALL SOLVE (J, FAIL, INV, JP, NNP) 
!
!   Upon failure issue message; take corrective action if possible
!
      IF (FAIL) THEN 
         IF (MYID == 0) WRITE (*, 300) NP(J), NH(J), METHOD(J) 
         IF (METHOD(J) /= 2) THEN 
            METHOD(J) = 2 
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            CALL ORTHSC 
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN 
!               CALL MATRIX 
!               CALL NEWCO (WTAEV) 
!            ENDIF 
            CALL SETLAG (EOL) 
            GO TO 1 
         ELSE 
            IF (MYID == 0) WRITE (*, 301) 
            !CALL TIMER (0)
            STOP  
         ENDIF 
      ENDIF 
!
!   Compute norm of radial function
!
      TA(1) = 0.D0 
      TA(2:MTP0) = (P(2:MTP0)**2+Q(2:MTP0)**2)*RP(2:MTP0) 
      MTP = MTP0 
 
      CALL QUAD (DNORM) 
 
!   Determine self-consistency [multiplied by SQRT(UCF(J))]
 
      CALL CONSIS (J) 
!
!   Normalize
!
      DNFAC = 1.D0/SQRT(DNORM) 
      P0 = P0*DNFAC 
      P(:MTP0) = P(:MTP0)*DNFAC 
      Q(:MTP0) = Q(:MTP0)*DNFAC 
!
!   Check if different method should be used or if improvement
!   count should be reduced
!
      DEL1 = ABS(1.D0 - ED2/E(J)) 
      IF (METHOD(J) == 1) THEN 
         DEL2 = MAX(ABS(1.D0 - SQRT(DNORM)),ABS(DNFAC - 1.D0)) 
         IF (DEL1<P005 .AND. DEL2>P2) THEN 
            METHOD(J) = 2 
            GO TO 1 
         ENDIF 
      ELSE 
         IF (DEL1<P0001 .AND. NSIC>1) NSIC = NSIC - 1 
      ENDIF 
!
!   Damp the orbital --- if not converged
!
      IF (SCNSTY(J) > ACCY) THEN 
         CALL DAMPCK (IPR, J, ED1, ED2) 
         ODAMPJ = ABS(ODAMP(J)) 
      ELSE 
         ODAMPJ = 0.D0                           ! take the whole new orbital 
      ENDIF 
      CALL DAMPOR (J, INV, ODAMPJ) 
 
!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.
 
      IF (ORTHST) THEN 
         !CALL orthor (J, inv)
         NWWW = NW 
         CALL ORTHY (NWWW, J, LSORT) 
      ENDIF 
!
!   Print details of iteration
!
      IF (MYID == 0)                                                &
         WRITE (*, 302) NP(J),NH(J),E(J),METHOD(J),PZ(J),SCNSTY(J), &
!cjb DNORM-1 -> SQRT(DNORM)-1
!cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP 
                        SQRT(DNORM)-1,ODAMPJ,JP,MF(J),INV,NNP
      DAMPMX = MAX(DAMPMX,ABS(ODAMPJ)) 
 
  300 FORMAT(/,' Failure; equation for orbital ',1I2,1A2,&
         ' could not be solved using method ',1I1) 
  301 FORMAT(/,/,' ****** Error in SUBROUTINE IMPROV ******'/,&
         ' Convergence not obtained'/) 
  302 FORMAT (1X,1I2,1A2,1P,1D16.7,1x,1I2,D11.3,1D10.2,1D10.2,&
              0P,F6.3,1x,1I5,1x,1I5,1x,1I2,1x,1I2)
      RETURN  
      END SUBROUTINE IMPROV 
