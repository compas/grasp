!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION GETYN ()
!                                                                      *
!   This  subprogram reads a response on  the default input unit; the  *
!   response must be either 'y' or 'n'. GETYN is .TRUE. if 'y' is en-  *
!   tered and .FALSE. if 'n' is entered.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 27 Aug 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:16   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE IOUNIT_C,     ONLY: ISTDI, ISTDE
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER (LEN = 1) :: RSPNS
!-----------------------------------------------
    1 CONTINUE
      READ (ISTDI, '(A)') RSPNS
      IF (RSPNS == 'y' .OR. RSPNS == 'Y') THEN
         GETYN = .TRUE.
      ELSE IF (RSPNS == 'n' .OR. RSPNS == 'N') THEN
         GETYN = .FALSE.
      ELSE
         WRITE (ISTDE, *) 'Expecting <y><cr> or <n><cr> ...'
         GO TO 1
      ENDIF
!
      RETURN
      END FUNCTION GETYN
