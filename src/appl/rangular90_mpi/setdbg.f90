!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBG (DEBUG, fullname)
!                                                                      *
!   This routine opens the  .dbg  file and sets the arrays that con-   *
!   trol debug printout from the GENMCP program.                       *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
!                                                                      *
!   Written by Farid A Parpia               Last update: 08 Dec 1992   *
!   Modified by Xinghong He                 Last update: 03 Jul 1998   *
!
!   Used by mcpvu, mcpmpivu
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:11:16  12/23/06
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE DEBUG_C
      USE DEFAULT_C, ONLY: NDEF
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(out) :: debug
      Character(LEN=*) :: fullname
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IERR
      LOGICAL :: YES
!
      CHARACTER(LEN = LEN (fullname)) FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
!
!   Initialise the arrays that control the debug printout
!
      DO I = 1, 5
         LDBPA(I) = .FALSE.
         LDBPG(I) = .FALSE.
      ENDDO
!
      RETURN      ! Skip these questions and return immediately
!
      IF (NDEF .EQ. 0) THEN
         RETURN
      ENDIF

      WRITE (istde,*) 'Generate debug printout?'
      DEBUG = GETYN ()
      IF (DEBUG) THEN
!
!   The  .dbg  file is formatted; open it on unit 99
!
         FORM = 'FORMATTED'
         STATUS = 'NEW'
!
         WRITE (istde,*) 'File ', fullname,' will be created as the'    &
     &,                 ' GENMCP DeBuG Printout File; '
         WRITE (istde,*) ' enter another file name if this is not '     &
     &,                 'acceptable; null otherwise:'
         READ (*,'(A)') FILNAM
!
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = fullname
!
    4    CALL OPENFL (99, FILNAM, FORM, STATUS, IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE (istde,*) 'Enter a name for the GENMCP DeBuG Printout'&
     &,                    ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
!
!   Set options for general printout
!
         WRITE (istde,*) ' Print out the machine constants used?'
         LDBPG(1) = GETYN ()
!
!   Set options for angular modules
!
         WRITE (istde,*) ' Printout from LODCSL? (Not used)'
         LDBPA(1) = GETYN ()
         WRITE (istde,*) ' Print out T coefficients?'
         LDBPA(2) = GETYN ()
         WRITE (istde,*) ' Print out Coulomb V coefficients?'
         LDBPA(3) = GETYN ()
         WRITE (istde,*) ' Print out sparse matrix definition arrays?'
         LDBPA(4) = GETYN ()

      ENDIF
!
      RETURN
      END SUBROUTINE SETDBG
