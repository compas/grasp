

!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBG
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
!                                                                      *
!   Written by Farid A Parpia               Last update: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE DEBUG_C
      USE DEFAULT_C
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IERR
      LOGICAL :: YES
      CHARACTER :: FILNAM*256, DEFNAM*11, FORM*11, STATUS*3
!-----------------------------------------------
!
!
!   Initialise the arrays that control the debug printout
!
      LDBPA = .FALSE.
!
      LDBPG = .FALSE.
!
      LDBPR = .FALSE.
!
      IF (NDEF == 0) RETURN

      WRITE (ISTDE, *) 'Generate debug printout?'
      YES = GETYN()
      IF (YES) THEN
!
!   The  .dbg  file is formatted; open it on unit 99
!
         DEFNAM = 'erwf.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
!
         WRITE (ISTDE, *) 'File  erwf.dbg  will be created as the ', &
            'ERWF DeBuG Printout File; '
         WRITE (ISTDE, *) 'enter another file name if this is not ', &
            'acceptable; null otherwise:'
         READ (*, '(A)') FILNAM
!
         IF (LEN_TRIM(FILNAM) == 0) FILNAM = DEFNAM
!
    4    CONTINUE
         CALL OPENFL (99, FILNAM, FORM, STATUS, IERR)
         IF (IERR /= 0) THEN
    5       CONTINUE
            WRITE (ISTDE, *) 'Enter a name for the ERWF DeBuG Printout ', &
               'file that is to be created:'
            READ (*, '(A)') FILNAM
            IF (LEN_TRIM(FILNAM) == 0) GO TO 5
            GO TO 4
         ENDIF
!
!   Set options for general printout
!
         WRITE (ISTDE, *) ' Print out the machine constants used?'
         YES = GETYN()
         IF (YES) LDBPG(1) = .TRUE.
         WRITE (ISTDE, *) ' Print out the physical constants used?'
         YES = GETYN()
         IF (YES) LDBPG(2) = .TRUE.
!
!   Set options for radial modules
!
         WRITE (ISTDE, *) ' Printout from RADGRD?'
         YES = GETYN()
         IF (YES) LDBPR(1) = .TRUE.
         WRITE (ISTDE, *) ' Printout from NUCPOT?'
         YES = GETYN()
         IF (YES) LDBPR(2) = .TRUE.
         WRITE (ISTDE, *) ' Printout from TFPOT?'
         YES = GETYN()
         IF (YES) LDBPR(26) = .TRUE.
!
!   Set options for angular modules
!
         WRITE (ISTDE, *) ' Printout from LODCSL?'
         YES = GETYN()
         IF (YES) LDBPA(1) = .TRUE.
!
      ENDIF
!
      RETURN
      END SUBROUTINE SETDBG
