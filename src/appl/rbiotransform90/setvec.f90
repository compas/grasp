!
!     ------------------------------------------------------------------
!     S E T V E C
!     ------------------------------------------------------------------
!
      SUBROUTINE SETVEC(VECTOR, VALUE, NDIM)
!
! VECTOR (*) = VALUE
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NDIM
      REAL(DOUBLE), INTENT(IN) :: VALUE
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(OUT) :: VECTOR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------
!
      VECTOR(:NDIM) = VALUE
!
      RETURN
      END SUBROUTINE SETVEC
