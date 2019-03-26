!***********************************************************************
!                                                                      *
      SUBROUTINE DEFCOR(J)
!                                                                      *
!   Compute the deferred corrections for orbital J .                   *
!                                                                      *
!                                          Last updated: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE grid_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: W3 = 1.0D00/120.0D00
      REAL(DOUBLE), PARAMETER :: W2 = -15.0D00*W3
      REAL(DOUBLE), PARAMETER :: W1 = 40.0D00*W3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, MFJM3, MFJM2
      LOGICAL :: FIRST
!-----------------------------------------------
!
      DATA FIRST/ .TRUE./
!
!   The deferred corrections for the first two points are
!   unnecessary, because the integration always commences
!   from the fourth point of the grid; this requires the
!   deferred correction at the third and subsequent points
!   only
!
      IF (FIRST) THEN
         DP(:2) = 0.0D00
         DQ(:2) = 0.0D00
         FIRST = .FALSE.
      ENDIF
!
!   Intermediate points
!
      MFJM3 = MF(J) - 3
      DO I = 3, MFJM3
!
         DP(I) = W3*(PF(I + 3,J) - PF(I - 2,J)) + &
                 W2*(PF(I + 2,J) - PF(I - 1,J)) + &
                 W1*(PF(I + 1,J) - PF(I,J))
!
         DQ(I) = W3*(QF(I + 3,J) - QF(I - 2,J)) + &
                 W2*(QF(I + 2,J) - QF(I - 1,J)) + &
                 W1*(QF(I + 1,J) - QF(I,J))
!
      END DO
!
!   Set remaining deferred corrections to zero: slopes are
!   small in this region
!
      MFJM2 = MF(J) - 2
      DP(MFJM2:N) = 0.0D00
      DQ(MFJM2:N) = 0.0D00
!
      RETURN
      END SUBROUTINE DEFCOR
