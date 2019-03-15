!***********************************************************************
!                                                                      *
      SUBROUTINE SETSUM(FULLNAME)
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]:  OPENFL.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Nov 1992   *
!   Modified by Xinghong He               Last revision:  3 Jul 1998   *
!
!  File shared by mcpblk, mcpmpi
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:54:22   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT(IN) :: FULLNAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR
      CHARACTER(LEN=120) :: FILNAM  ! The original did not compile
      CHARACTER(LEN=11) :: FORM
      CHARACTER(LEN=3)  :: STATUS
!-----------------------------------------------

      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      WRITE (ISTDE, *) 'File ', FULLNAME, ' will be created as the', &
         ' GENMCP SUMmary File;'
      WRITE (ISTDE, *) 'enter another file name if this is not ', &
         'acceptable; null otherwise:'
      READ (*, '(A)') FILNAM
!
      IF (LEN_TRIM(FILNAM) == 0) FILNAM = FULLNAME
!
    1 CONTINUE
      CALL OPENFL (24, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
    2    CONTINUE
         WRITE (ISTDE, *) 'Enter a name for the GENMCP SUMmary', &
            ' File that is to be created:'
         READ (*, '(A)') FILNAM
         IF (LEN_TRIM(FILNAM) == 0) GO TO 2
         GO TO 1
      ENDIF
!
      RETURN
      END SUBROUTINE SETSUM
