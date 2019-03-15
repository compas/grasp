!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION HOVLAP (P, Q, MTPO, NP, KAPPA, Z)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
      INTEGER, INTENT(IN) :: MTPO
      INTEGER  :: NP
      INTEGER  :: KAPPA
      REAL(DOUBLE)  :: Z
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: P
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: Q
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTPH, I
      REAL(DOUBLE) :: EH, PZH, RESULT
!-----------------------------------------------
!
!   Set up the hydrogenic orbital
!
      CALL DCBSRW (NP, KAPPA, Z, EH, PZH, PH, QH, MTPH)
!
!   Compute the overlap
!
      MTP = MIN(MTPH,MTPO)
      TA(1) = 0.0D00
      TA(2:MTP) = (P(2:MTP)*PH(2:MTP)+Q(2:MTP)*QH(2:MTP))*RP(2:MTP)
      CALL QUAD (RESULT)
!
      HOVLAP = RESULT
!
      RETURN
      END FUNCTION HOVLAP
