!
!
!     ------------------------------------------------------------------
!     S C A L V E
!     ------------------------------------------------------------------
!
      SUBROUTINE SCALVE(VECTOR, FACTOR, NDIM) 
!
! CALCULATE SCALAR(FACTOR) TIMES VECTOR
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
      REAL(DOUBLE), INTENT(IN) :: FACTOR 
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(INOUT) :: VECTOR
!-----------------------------------------------
!
      VECTOR(:NDIM) = VECTOR(:NDIM)*FACTOR 
!
      RETURN  
      END SUBROUTINE SCALVE 
