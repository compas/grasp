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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
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
      CHARACTER :: FILNAM1*256, FILNAM2*256, DEFNAM*11, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File  hfs92.sum  is FORMATTED
!
      K = INDEX(NAME,' ')
      IF (NCI == 0) THEN
         FILNAM1 = NAME(1:K-1)//'.ch'
         FILNAM2 = NAME(1:K-1)//'.choffd'
      ELSE
         FILNAM1 = NAME(1:K-1)//'.h'
         FILNAM2 = NAME(1:K-1)//'.hoffd'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      CALL OPENFL (29, FILNAM1, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM1
         STOP
      ENDIF
!
      CALL OPENFL (24, FILNAM2, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM2
         STOP
      ENDIF
!
      RETURN
      END SUBROUTINE SETSUM
