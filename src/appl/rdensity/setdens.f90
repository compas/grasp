!***********************************************************************
!                                                                      *
      SUBROUTINE SETDENS(NAME, NCI) 
!                                                                      *
!   Open the  .d  files on stream 35                                   *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by J. Ekman                                  23 Nov 2013   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NCI 
      CHARACTER, INTENT(IN) :: NAME*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: K, IERR 
      CHARACTER :: FILNAM*256, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File  rdensity.sum  is FORMATTED
!
      K = INDEX(NAME,' ') 
      IF (NCI == 0) THEN 
         FILNAM = NAME(1:K-1)//'.cd' 
      ELSE 
         FILNAM = NAME(1:K-1)//'.d' 
      ENDIF 
      FORM = 'FORMATTED' 
      STATUS = 'NEW' 
!
      CALL OPENFL (35, FILNAM, FORM, STATUS, IERR) 
      IF (IERR /= 0) THEN 
         WRITE (6, *) 'Error when opening', FILNAM
         STOP  
      ENDIF 
!
      RETURN  
      END SUBROUTINE SETDENS
