      SUBROUTINE DINIT(N, A, X, INCX) 
!     ==================================================================
!
!     PURPOSE ... INITIALIZES REAL*8           VECTOR TO
!                 A CONSTANT VALUE 'A'
!
!     CREATED ... APR. 14, 1987
!
!     ==================================================================
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:18   2/14/04  
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
      INTEGER, INTENT(IN) :: N 
      INTEGER, INTENT(IN) :: INCX 
      REAL(DOUBLE), INTENT(IN) :: A 
      REAL(DOUBLE), DIMENSION(*), INTENT(OUT) :: X
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: XADDR, I 
!-----------------------------------------------
!
      IF (INCX == 1) THEN 
!
!         ----------------------------------
!         ... UNIT INCREMENT (STANDARD CASE)
!         ----------------------------------
!
         X(:N) = A 
!
      ELSE 
!
!         ----------------------
!         ... NON-UNIT INCREMENT
!         ----------------------
!
         XADDR = 1 
         IF (INCX < 0) XADDR = ((-N) + 1)*INCX + 1 
!
         X(XADDR:(N-1)*INCX+XADDR:INCX) = A 
!
      ENDIF 
!
      RETURN  
!
      END SUBROUTINE DINIT 
