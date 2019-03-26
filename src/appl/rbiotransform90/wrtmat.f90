!***********************************************************************
!                                                                      *
      SUBROUTINE WRTMAT(A, NROW, NCOL, NMROW, NMCOL)
!                                                                      *
!***********************************************************************
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
      INTEGER, INTENT(IN) :: NROW
      INTEGER, INTENT(IN) :: NCOL
      INTEGER, INTENT(IN) :: NMROW
      INTEGER, INTENT(IN) :: NMCOL
      REAL(DOUBLE), INTENT(IN) :: A(NMROW,NMCOL)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J
!-----------------------------------------------

      DO I = 1, NROW
         WRITE (6, 1010) I, (A(I,J),J=1,NCOL)
      END DO

 1010 FORMAT('0',I5,2X,4(1X,E14.7),/,(' ',7X,4(1X,E14.7)))

      RETURN
      END SUBROUTINE WRTMAT
