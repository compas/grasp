!
!     ------------------------------------------------------------------
!     T R P M A T
!     ------------------------------------------------------------------
!
      SUBROUTINE TRPMAT(XIN, NROW, NCOL, XOUT) 
!
! XOUT(I,J) = XIN(J,I)
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
      INTEGER , INTENT(IN) :: NROW 
      INTEGER , INTENT(IN) :: NCOL 
      REAL(DOUBLE) , INTENT(IN) :: XIN(NROW,NCOL) 
      REAL(DOUBLE) , INTENT(OUT) :: XOUT(NCOL,NROW) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IROW, ICOL 
!-----------------------------------------------
!
      XOUT = TRANSPOSE(XIN) 
!
      RETURN  
      END SUBROUTINE TRPMAT 
