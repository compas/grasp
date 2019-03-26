!***********************************************************************
!                                                                      *
      SUBROUTINE OPENFL(NFILE, FILNAM, RFORM, RSTAT, IERR)
!                                                                      *
!   Issues OPEN for file with unit number NFILE, name FILNAM, format   *
!   RFORM, status RSTAT.  If this is successful the head is position-  *
!   ed to the beginning of the file and IERR is 0; otherwise IERR is   *
!   set to 1.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]:  none.                                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:04   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE IOUNIT_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NFILE
      INTEGER, INTENT(OUT) :: IERR
      CHARACTER (LEN = *), INTENT(IN) :: FILNAM
      CHARACTER (LEN = *), INTENT(IN) :: RFORM
      CHARACTER (LEN = *), INTENT(IN) :: RSTAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOS, LOC
!-----------------------------------------------
!
!
      OPEN(NFILE, FILE=FILNAM, FORM=RFORM, STATUS='UNKNOWN', IOSTAT=IOS, &
         POSITION='asis')
!
      IF (IOS == 0) THEN
         REWIND (NFILE)
         IERR = 0
      ELSE
         LOC = LEN_TRIM(FILNAM)
         WRITE (ISTDE, *) 'OPENFL: Error opening file ', FILNAM(1:LOC), ' as '&
            , RSTAT, ';'
         WRITE (ISTDE, *) 'The argument RSTAT=', RSTAT, ' is not used !'
         IERR = 1
      ENDIF
!
      RETURN
      END SUBROUTINE OPENFL
