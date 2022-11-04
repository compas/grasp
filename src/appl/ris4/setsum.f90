!***********************************************************************
!                                                                      *
      SUBROUTINE SETSUM(NAME, NCI)
!                                                                      *
!   Open the  .sum  files on stream 24 and 29                          *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!                                                                      *
!   Updated by Per Jonsson                               28 Oct 1999   *
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
         FILNAM = NAME(1:K-1)//'.ci'
      ELSE
         FILNAM = NAME(1:K-1)//'.i'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      CALL OPENFL (24, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM
         STOP
      ENDIF
!
      RETURN
      END SUBROUTINE SETSUM
