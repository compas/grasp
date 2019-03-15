!***********************************************************************
!                                                                      *
      SUBROUTINE QED(JSTATE, SLFINT, UCF)
!                                                                      *
!   This  routine estimates corrections to the  energy levels due to   *
!   self-energy.                                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SCREEN.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!                                           Last update: 30 Oct 1992   *
!   Modified by Xinghong He                 Last update: 24 Jun 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW, NNNP
      USE def_C
      USE eigv_C
      !USE grid_C
      !USE horb_C
      USE npar_C
      USE orb_C, ONLY: nw, ncf, np, nak
      USE prnt_C
      !USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      USE ratden_I
      USE fzalf_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JSTATE
      REAL(DOUBLE), INTENT(OUT) :: SLFINT(NNNW)
      REAL(DOUBLE), INTENT(OUT) :: UCF(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MAXITER, J, I, II, NPJ, KAPPA, MFJ
      REAL(DOUBLE), DIMENSION(NNNP) :: PTEMP, QTEMP
      REAL(DOUBLE) :: UCFJ, ZEFF, RATIO, VALU
      CHARACTER :: NPCHAR, NAKCHAR*2
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      MAXITER = 20
!
! Modified so that UCFJ describes the current eigenstate
!
      DO J = 1, NW
         UCFJ = 0.0D00
!         DO 3 I = 1,NVEC
         I = JSTATE
         DO II = 1, NCF
            UCFJ = UCFJ + DBLE(IQ(J,II))*EVEC(II + (I - 1)*NCF)**2
         END DO
!    3    CONTINUE
!     print *, ucfj,'ucf'
         UCF(J) = UCFJ
! zou    UCF(J) = UCFJ/DBLE (NCF)
      END DO
!
      DO J = 1, NW
!
         NPJ = NP(J)
!
         IF (NPJ <= 8) THEN
!
!   Only orbitals with principal quantum number 8 or less can
!   be treated by this section of code
!
            KAPPA = NAK(J)
!
!   Begin by transferring the function to a temporary array
!
            MFJ = MF(J)
!
            PTEMP(1) = 0.0D00
            QTEMP(1) = 0.0D00
            DO I = 2, MFJ
               PTEMP(I) = PF(I,J)
               QTEMP(I) = QF(I,J)
            END DO
!
            ZEFF = Z
            RATIO = RATDEN(PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZEFF)
            VALU = RATIO*FZALF(NPJ,KAPPA,ZEFF)/DBLE(NPJ**3)
!
            SLFINT(J) = VALU*ZEFF**4/(PI*C**3)
!        print *, 'No. orb.=',j,' Zeff = ',zeff
!    &   , 'Scale= ',ratio
!    &   , 'S.E. = ',slfint(j)*2*13.6058,slfint(j)/ratio*2*13.6058
!
         ELSE
!
!   The self-energy for orbitals with principal quantum number
!   greater than 8 is set to zero
!
            SLFINT(J) = 0.0D00
!
         ENDIF
!
      END DO
!
!   Deallocate storage for the `generalised occupation numbers'
!
!
      RETURN
      END SUBROUTINE QED


!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RATDEN (P, Q, MTPO, NP, KAPPA, Z)
!                                                                      *
!   This subprogram computes the overlap of the orbital tabulated in   *
!   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
!   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Switches:
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE grid_C
      USE horb_C, ONLY: ph, qh
      USE tatb_C, ONLY: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dcbsrw_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MTPO
      INTEGER  :: NP
      INTEGER  :: KAPPA
      REAL(DOUBLE)  :: Z
      REAL(DOUBLE), INTENT(IN) :: P(NNNP)
      REAL(DOUBLE), INTENT(IN) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTPH, I, K
      REAL(DOUBLE) :: EH, PZH, RESULT, RESULT1
!-----------------------------------------------
!
!
!
!   Set up the hydrogenic orbital
!
      CALL DCBSRW (NP, KAPPA, Z, EH, PZH, PH, QH, MTPH)
!
!   Compute the overlap
!
      MTP = MIN(MTPH,MTPO)
      DO I = 2, MTP
         IF (RP(I) > 0.0219) CYCLE
         K = I
      END DO
      MTP = K
      TA(1) = 0.0D00
      TA(2:MTP) = (P(2:MTP)*P(2:MTP)+Q(2:MTP)*Q(2:MTP))*RP(2:MTP)
!         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
      CALL QUAD (RESULT)
      TA(2:MTP) = (PH(2:MTP)*PH(2:MTP)+QH(2:MTP)*QH(2:MTP))*RP(2:MTP)
      CALL QUAD (RESULT1)
!
      RATDEN = RESULT/RESULT1
!
      RETURN
      END FUNCTION RATDEN
