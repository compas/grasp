!***********************************************************************     
      SUBROUTINE OUTSDA(LPRINT,NNONZ, NCF) 
!***********************************************************************     
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      IMPLICIT NONE 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NNONZ 
      INTEGER, INTENT(IN) :: NCF 
      LOGICAL, INTENT(IN) :: LPRINT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=9), PARAMETER :: MYNAME ='outsdampi'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NNONZ_A 
!-----------------------------------------------
      NNONZ_A = NNONZ 
 
      WRITE(6, *)' ... complete; density of non-zero elements of H(DC):', &
         NNONZ_A, '/', (NCF*(NCF + 1))/2 
!
!   Debug printout
!
      IF (LPRINT) THEN 
         WRITE (99, *) 
         WRITE (99, *) 'From ', MYNAME, ' :' 
         WRITE (99, 301) NNONZ_A 
 
         WRITE (6, *) 'This part not finished. See ', MYNAME 
         WRITE (99, *) 'This part not finished. See ', MYNAME 
      ENDIF 
 
  301 FORMAT(' Number of nonzero elements in H(DC): ',1I4) 
  302 FORMAT(' Column ',1I2,', row ',1I2,', sparse matrix index ',1I4) 
 
      RETURN  
      END SUBROUTINE OUTSDA 
