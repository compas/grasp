!***********************************************************************
!                                                                      *
      SUBROUTINE SETISO(FNAME)
!                                                                      *
!   Open, check, load data from and close the  .iso  file. This file   *
!   is always attached to stream  22 in  RSCF92,  RCI92,  HFS92, and   *
!   OSCL92.                                                            *
!   Filename - fname is moved to the argument list. This subroutine
!   is ready to replace the one in lib/lib92, but need to check
!   other application programs before doing do.
!                                                                      *
!   Call(s) to: [LIB92]: LODISO, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified    by Xinghong He            Last revision:  1 Jun 1998   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:35   2/14/04
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
      USE lodiso_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = *), INTENT(IN) :: FNAME
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=3), PARAMETER  :: STATUS = 'OLD'
      CHARACTER(LEN=6), PARAMETER  :: MYNAME = 'SETISO'
      CHARACTER(LEN=9), PARAMETER  :: FORM = 'FORMATTED'
      CHARACTER(LEN=14), PARAMETER :: SIGNATURE = 'Atomic number:'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER           :: LENF, IERR, IOS
      LOGICAL           :: FOUND
      CHARACTER(LEN=14) :: STR
!-----------------------------------------------
!
      INQUIRE(FILE=FNAME, EXIST=FOUND)
      IF (.NOT.FOUND) THEN
         LENF = LEN_TRIM(FNAME)
         WRITE (ISTDE,*)                                              &
                     MYNAME,'- file: ',FNAME(1:LENF),' does not exist'
         STOP
      ENDIF
!
      CALL OPENFL (22, FNAME, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) 'Error opening isodata file: ', FNAME(1:LENF)
         STOP
      ENDIF
!
!   Check the first record of the file; if not as expected, try again
!
      READ (22, '(A)', IOSTAT=IOS) STR
      IF (IOS/=0 .OR. STR/=SIGNATURE) THEN
         WRITE (ISTDE, *) 'Not an ISOtope Data File;'
         CLOSE(22)
         STOP
      ENDIF
!
!   Load data from the .iso file and then close it.
!
      CALL LODISO
      CLOSE(22)

      RETURN
      END SUBROUTINE SETISO
