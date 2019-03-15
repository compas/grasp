!***********************************************************************
!                                                                      *
      SUBROUTINE SETCSLA(NAME, NCORE, IGG)
!-----------------------------------------------
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LENGTH, LODCSL, OPENFL.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodcslBio_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE, IGG
      CHARACTER, INTENT(IN) :: NAME*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IERR, IOS
      LOGICAL :: FOUND
      CHARACTER :: FILNAM*256, RECORD*15, DEFNAM*11, FORM*11, STATUS*3
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      K = INDEX(NAME,' ')
      FILNAM = NAME(1:K-1)//'.c'
      FORM = 'FORMATTED'
      STATUS = 'OLD'

!
      CALL OPENFL (21, FILNAM, FORM, STATUS, IERR)
      IF (IERR == 1) THEN
         WRITE (6, *) 'Error when opening', FILNAM
         STOP
      ENDIF
!
!   Check the first record of the file; if not as expected, try again
!
      READ (21, '(1A15)', IOSTAT=IOS) RECORD
      IF (IOS/=0 .OR. RECORD(1:15)/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         CLOSE(21)
      ENDIF
!
!   Load data from the  .csl  file
!
      CALL LODCSLBio (NCORE,IGG)
!
!   Close the  .csl  file
!
      CLOSE(21)
!
      RETURN
      END SUBROUTINE SETCSLA
