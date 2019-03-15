!     ------------------------------------------------------------------
!     P R S Y M
!     ------------------------------------------------------------------
!
      SUBROUTINE PRSYM(A, MATDIM)
! PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
! THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
!-----------------------------------------------
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
      INTEGER, INTENT(IN) :: MATDIM
      REAL(DOUBLE), DIMENSION(*), INTENT(IN) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JSTART, JSTOP, I, J
!-----------------------------------------------
      JSTART = 1
      JSTOP = 0
      DO I = 1, MATDIM
         JSTART = JSTART + I - 1
         JSTOP = JSTOP + I
         WRITE (6, 1010) I, (A(J),J=JSTART,JSTOP)
      END DO
      RETURN
 1010 FORMAT('0',2X,I3,5(1X,E14.7),/,(' ',5X,5(1X,E14.7)))
      RETURN
      END SUBROUTINE PRSYM
