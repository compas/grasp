

!***********************************************************************
!                                                                      *
      SUBROUTINE PRTREM(ALL)
!                                                                      *
!   Prints a list of subshells that remain to be estimated.  ALL  is   *
!   .TRUE.  if all subshells have been estimated and  .FALSE. other-   *
!   wise.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE LEFT_C
      USE ORB_C, ONLY: NW,  NP, NAK, NH
      USE IOUNIT_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL , INTENT(OUT) :: ALL
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /ORB2/
!     COMMON /ORB2/ NCF, NW, PNTRIQ
!     REAL   PNTRIQ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IEND, IBEG, LENTH
      CHARACTER(LEN=12) :: CNUM
      CHARACTER(LEN=80) :: RECORD*80
!-----------------------------------------------
!
!     Print *, ' Entering PRTREM'
!   Determine if there are any subshell radial wavefunctions that
!   remain to be estimated
!
      DO I = 1, NW
         IF (SET(I)) CYCLE
         ALL = .FALSE.
         GO TO 2
      END DO
      ALL = .TRUE.
      GO TO 4
!
!   Print a list of subshell radial wavefunctions that remain to
!   be estimated; this list is no more than 80 characters wide
!
    2 CONTINUE
      WRITE (ISTDE, *) 'The following subshell radial wavefunctions ', &
         'remain to be estimated:'
      IEND = 0
      DO I = 1, NW
         IF (SET(I)) CYCLE
         IF (IEND > 75) THEN
            WRITE (ISTDE, *) RECORD(1:IEND)
            IEND = 0
         ENDIF
         IF (IEND > 0) THEN
            IBEG = IEND + 1
            IEND = IBEG
            RECORD(IBEG:IEND) = ' '
         ENDIF
         IBEG = IEND + 1
         CALL CONVRT (NP(I), CNUM, LENTH)
         IEND = IBEG + LENTH - 1
         RECORD(IBEG:IEND) = CNUM(1:LENTH)
         IF (NAK(I) < 0) THEN
            LENTH = 1
         ELSE
            LENTH = 2
         ENDIF
         IBEG = IEND + 1
         IEND = IBEG + LENTH - 1
         RECORD(IBEG:IEND) = NH(I)(1:LENTH)
      END DO
      IF (IEND > 1) WRITE (ISTDE, *) RECORD(1:IEND)
!
    4 CONTINUE
      RETURN
      END SUBROUTINE PRTREM
