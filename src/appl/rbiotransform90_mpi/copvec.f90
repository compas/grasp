!
!     ------------------------------------------------------------------
!     C O P V E C
!     ------------------------------------------------------------------
!
      SUBROUTINE COPVEC(FROM, TO, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NDIM 
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(IN)  :: FROM
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(OUT) :: TO
!
      TO(:NDIM) = FROM(:NDIM) 
!
      RETURN  
      END SUBROUTINE COPVEC 
