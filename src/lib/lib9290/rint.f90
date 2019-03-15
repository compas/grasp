!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINT (I, J, K)
!                                                                      *
!   The value of RINT is an approximation to:                          *
!                                                                      *
!              k                                                       *
!         I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)       *
!                      i     j       i     j                           *
!                                                                      *
!   where   I ( G(r) ; range )  denotes  the  integral  of G(r) over   *
!   range.                                                             *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:22   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE GRID_C
      USE TATB_C, ONLY: TA, MTP
      USE WAVE_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I
      INTEGER  :: J
      INTEGER, INTENT(IN) :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------
!
!
!   Determine the maximum tabulation point for the integrand
!
      MTP = MIN(MF(I),MF(J))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
      TA(1) = 0.0D00
      DO L = 2, MTP
         TA(L) = R(L)**K*(PF(L,I)*PF(L,J) + QF(L,I)*QF(L,J))*RP(L)
      END DO
!
!   Perform the quadrature
!
      CALL QUAD (RESULT)
      RINT = RESULT
!
      RETURN
!
      END FUNCTION RINT
