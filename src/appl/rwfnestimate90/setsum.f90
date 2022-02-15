

!***********************************************************************
!                                                                      *
      SUBROUTINE SETSUM
      USE IOUNIT_C
!                                                                      *
!   Open the  .sum  file on stream 24.                                 *
!                                                                      *
!   Call(s) to: [LIB92]:  OPENFL.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR
      CHARACTER :: FILNAM*256, DEFNAM*11, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File  erwf.sum  is FORMATTED
!
      DEFNAM = 'erwf.sum'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      WRITE (ISTDE, *) 'File  erwf.sum  will be created as the ', &
         'ERWF SUMmary File; '
      WRITE (ISTDE, *) 'enter another file name if this is not ', &
         'acceptable; null otherwise:'
      READ (*, '(A)') FILNAM
!
      IF (LEN_TRIM(FILNAM) == 0) FILNAM = DEFNAM
!
    1 CONTINUE
      CALL OPENFL (24, FILNAM, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
    2    CONTINUE
         WRITE (ISTDE, *) 'Enter a name for the ERWF SUMmary File ', &
            'that is to be created:'
         READ (*, '(A)') FILNAM
         IF (LEN_TRIM(FILNAM) == 0) GO TO 2
         GO TO 1
      ENDIF
!
      RETURN
      END SUBROUTINE SETSUM
