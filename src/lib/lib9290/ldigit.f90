!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION LDIGIT (CST) 
!                                                                      *
!   .TRUE.  if  CST  is the ASCII representation of a decimal digit;   *
!   .FALSE. otherwise.                                                 *
!                                                                      *
!   Written by Farid A. Parpia             Last revised: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:52   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER , INTENT(IN) :: CST 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      CHARACTER, DIMENSION(0:9) :: CDGT 
!-----------------------------------------------
!
!
!
      DATA CDGT/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/  
!
      DO I = 0, 9 
         IF (CST /= CDGT(I)) CYCLE  
         LDIGIT = .TRUE. 
         GO TO 2 
      END DO 
      LDIGIT = .FALSE. 
!
    2 CONTINUE 
      RETURN  
      END FUNCTION LDIGIT 
