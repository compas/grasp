!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ITRIG (I1, I2, I3)
!                                                                      *
!   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
!   form a triangle; 0, otherwise.                                     *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:47   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I1
      INTEGER, INTENT(IN) :: I2
      INTEGER, INTENT(IN) :: I3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I4
!-----------------------------------------------
!
      I4 = I2 - I3
      IF (I1>=ABS(I4) + 1 .AND. I1<=I2+I3-1) THEN
         ITRIG = 1
      ELSE
         ITRIG = 0
      ENDIF
!
      RETURN
      END FUNCTION ITRIG
