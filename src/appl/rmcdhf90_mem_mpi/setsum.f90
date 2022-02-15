!***********************************************************************
      SUBROUTINE SETSUM(FILNAM)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
     USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: FILNAM*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*9, PARAMETER :: FORM = 'FORMATTED'
      CHARACTER*3, PARAMETER :: STATUS = 'NEW'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR
!-----------------------------------------------

      CALL OPENFL (24, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) 'Error when opening ', FILNAM
         STOP
      ENDIF

      RETURN
      END SUBROUTINE SETSUM
