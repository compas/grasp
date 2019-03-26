!***********************************************************************
!                                                                      *
      SUBROUTINE SETSUM(NAME, NCI)
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
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
      CHARACTER :: FILNAM*256, DEFNAM*11, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File  sms92.sum  is FORMATTED
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
