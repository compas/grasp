!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBG(DBGFILE)
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 10 Dec 1992   *
!   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:22:29   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE DEBUG_C
      USE default_C
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
      CHARACTER(LEN=*) , INTENT(IN) :: DBGFILE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*9, PARAMETER :: FORM = 'FORMATTED'
      CHARACTER*3, PARAMETER :: STATUS = 'NEW'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, LENDBG, LENFIL, IERR
!GG      CHARACTER(LEN=120) :: FILNAM
      CHARACTER(LEN = LEN (dbgfile)) :: FILNAM
!-----------------------------------------------
!
!   Initialise the arrays that control the debug printout
!   These serve as the default settings.
!
      LDBPA = .FALSE.
      LDBPG = .FALSE.

      LDBPR = .FALSE.
!
      IF (NDEF == 0) RETURN

!   Even in non-default, the user can choose not to have debug
!   print-out

      WRITE (ISTDE, '(A)') 'Generate debug output?  (y/n) '
      IF (.NOT.GETYN()) RETURN

      LENDBG = LEN_TRIM(DBGFILE)

      WRITE (ISTDE, *) 'File  ', DBGFILE(1:LENDBG), &
      '  will be                           created as the RSCF92 DeBuG Printout&
      & File;'
      WRITE (ISTDE, *) 'enter another file name if this is not', &
         ' acceptable; null otherwise:'

  123 CONTINUE
      READ (*, '(A)') FILNAM

      FILNAM = ADJUSTL(FILNAM)
      LENFIL = LEN_TRIM(FILNAM)
      IF (LENFIL == 0) THEN
         FILNAM = DBGFILE
      ELSE IF (LENFIL > LENDBG) THEN
         WRITE (ISTDE,*) 'File name too long, (> ', LENDBG, '); redo...'
         GO TO 123
      ENDIF

      CALL OPENFL (99, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) 'File name not accepted; redo...'
         GO TO 123
      ENDIF
!
!   Set options for general printout
!
      WRITE (ISTDE, *) 'Print out the machine constants used?'
      LDBPG(1) = GETYN()
      WRITE (ISTDE, *) 'Print out the physical constants used?'
      LDBPG(2) = GETYN()
      WRITE (ISTDE, *) 'Printout from FNDBLK?'
      LDBPG(3) = GETYN()
      WRITE (ISTDE, *) 'Print out the Hamiltonian matrix?'
      LDBPG(4) = GETYN()
      WRITE (ISTDE, *) 'Print out the eigenvectors?'
      LDBPG(5) = GETYN()
!      LDBPG(1:5) = .TRUE.
!
!   Set options for printout from radial modules
!
      WRITE (ISTDE, *) 'Printout from RADGRD?'
      LDBPR(1) = GETYN()
      WRITE (ISTDE, *) 'Printout from NUCPOT?'
      LDBPR(2) = GETYN()
      WRITE (ISTDE, *) 'Printout from LODRWF?'
      LDBPR(3) = GETYN()
      WRITE (ISTDE, *) 'Print out I(ab) integrals?'
      LDBPR(4) = GETYN()
      WRITE (ISTDE, *) 'Print out Slater integrals?'
      LDBPR(10) = GETYN()
      WRITE (ISTDE, *) 'Make summary printout on progress', &
         ' of each iteration in SOLVE?'
      LDBPR(22) = GETYN()
      WRITE (ISTDE, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions on', ' each iteration in SOLVE?'
      LDBPR(23) = GETYN()
      WRITE (ISTDE, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions', ' after each SCF cycle?'
      LDBPR(24) = GETYN()
      WRITE (ISTDE, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions on', ' convergence?'
      LDBPR(25) = GETYN()
      WRITE (ISTDE, *) 'List compositions of exchange', ' potentials?'
      LDBPR(27) = GETYN()
      WRITE (ISTDE, *) 'Tabulate and make printer plots', &
         ' of exchange potentials?'
      LDBPR(28) = GETYN()
      WRITE (ISTDE, *) 'List compositions of direct', ' potentials?'
      LDBPR(29) = GETYN()
      WRITE (ISTDE, *) 'Tabulate and make printer plots', &
         ' of direct potentials?'
      LDBPR(30) = GETYN()
!      LDBPR(1:30) = .TRUE.
!
!   Set options for printout of angular coefficients
!
      WRITE (ISTDE, *) ' Printout from LODCSL?'
      LDBPA(1) = GETYN()
      WRITE (ISTDE, *) ' Print out T coefficients?'
      LDBPA(2) = GETYN()
      WRITE (ISTDE, *) ' Print out V coefficients?'
      LDBPA(3) = GETYN()
!      LDBPA(1:3) = .TRUE.

      RETURN
      END SUBROUTINE SETDBG
