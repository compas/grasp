!***********************************************************************
!                                                                      *
      SUBROUTINE QED_SLFEN(SLFINT) 
!                                                                      *
!   This  routine estimates the F(Z\alpha) function of self energy for *
!   each orbital.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, IQ, QUAD, RALLOC, SHIELD.      *
!               [RCI92]: FZALF, HOVLAP.                                *
!                                                                      *
!   Modified from subroutine QED by Yu Zou, Last update: 13 Mar 2000   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE parameter_def,   ONLY: NNNW, NNNP
      USE def_C
      USE eigv_C
      USE npar_C
      USE orb_C, ONLY: nw, np, nak
      USE prnt_C
      USE wave_C
      USE qedcut_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ratden_I 
      USE fzalf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW), INTENT(OUT) :: SLFINT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MAXITER, J, NPJ, KAPPA, MFJ, I, NPJMAX
!GG      REAL(DOUBLE), DIMENSION(1) :: UCF 
      REAL(DOUBLE), DIMENSION(NNNP) :: PTEMP, QTEMP 
      REAL(DOUBLE) :: ZEFF, RATIO, VALU 
      CHARACTER :: NPCHAR, NAKCHAR*2 
!-----------------------------------------------
!
! Pre-set tolerable number for iteration in finding effective
! nuclear charge.
!
      MAXITER = 20 
!Per
      IF (NQEDCUT.EQ.1) THEN
         NPJMAX = NQEDMAX
      ELSE
         NPJMAX = 8
      END IF
!Per
!
      DO J = 1, NW 
!
         NPJ = NP(J) 
!
         IF (NPJ <= NPJMAX) THEN 
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
            ZEFF = Z 
            RATIO = RATDEN(PTEMP,QTEMP,MFJ,NPJ,KAPPA,ZEFF) 
            VALU = RATIO*FZALF(NPJ,KAPPA,ZEFF)/DBLE(NPJ**3) 
            SLFINT(J) = VALU*ZEFF**4/(PI*C**3) 
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
      END SUBROUTINE QED_SLFEN 
