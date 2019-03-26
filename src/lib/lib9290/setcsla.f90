!***********************************************************************
!                                                                      *
      SUBROUTINE SETCSLA(NAME, NCORE)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LODCSL, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:33   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodcsl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      CHARACTER (LEN = 24), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: K, IERR, IOS
      CHARACTER(LEN=3)   :: STATUS
      CHARACTER(LEN=9)   :: FORM
      CHARACTER(LEN=15)  :: RECORD
      CHARACTER(LEN=256) :: FILNAM
!-----------------------------------------------
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
         WRITE (ISTDE, *) 'Error when opening', FILNAM
         STOP
      ENDIF
!
!   Check the first record of the file; if not as expected, try again
!
      READ (21, '(1A15)', IOSTAT=IOS) RECORD
      IF (IOS/=0 .OR. RECORD(1:15)/='Core subshells:') THEN
         WRITE (ISTDE, *) 'Not a Configuration Symmetry List File;'
         CLOSE(21)
      ENDIF
!
!   Load data from the  .csl  file
!
      CALL LODCSL (NCORE)
!
!   Close the  .csl  file
!
      CLOSE(21)
!
      RETURN
      END SUBROUTINE SETCSLA
