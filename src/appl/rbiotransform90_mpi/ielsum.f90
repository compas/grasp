!
!     ------------------------------------------------------------------
!     I E L S U M
!     ------------------------------------------------------------------
!
      INTEGER FUNCTION IELSUM (IVEC, NELMNT) 
!
! Sum elements of integer array
!
!************************************************************************     
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NELMNT 
      INTEGER, INTENT(IN) :: IVEC(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISUM, IEL 
!-----------------------------------------------
!
      ISUM = 0 
      ISUM = SUM(IVEC(:NELMNT)) 
!
      IELSUM = ISUM 
!
      RETURN  
      END FUNCTION IELSUM 
