!***********************************************************************
!                                                                      *
      SUBROUTINE PRSRCN(RECORD, NCORE, IOCCS, IERR)
!                                                                      *
!   READs and parses a string that specifies a configuration.          *
!                                                                      *
!   Written by Farid A Parpia              Last revised: 16 Oct 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:15   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE ORB_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)           :: NCORE
      INTEGER, INTENT(OUT)          :: IERR
      CHARACTER(LEN=256), INTENT(IN) :: RECORD
      INTEGER, DIMENSION(NNNW), INTENT(OUT) :: IOCCS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,ISTART,IEND,LENTH,IOS,NPI,J,ISHELL,IOSTRT,IOEND, &
         IOCCI,NKJI
      CHARACTER(LEN=2) :: SYMI
      CHARACTER(LEN=5) :: FORM
      CHARACTER(LEN=1), DIMENSION(3) :: CNUM = (/'1', '2', '3'/)
      CHARACTER :: RECI
!-----------------------------------------------
!
!   Initialise IOCCS for the peel subshells
!
      IOCCS(NCORE+1:NW) = 0
!
!   Parse RECORD from left to right
!
      ISTART = 0
      I = 1
    2 CONTINUE
      RECI = RECORD(I:I)
      IF (RECI == '(') THEN
         IEND = I - 1
    3    CONTINUE
         RECI = RECORD(IEND:IEND)
         IF (RECI == ' ') THEN
            IEND = IEND - 1
            GO TO 3
         ELSE IF (RECI == '-') THEN
            READ (RECORD(IEND-1:IEND), '(1A2)') SYMI
            IEND = IEND - 2
         ELSE
            SYMI(2:2) = ' '
            READ (RECI, '(1A1)') SYMI(1:1)
            IEND = IEND - 1
         ENDIF
         LENTH = IEND - ISTART + 1
         FORM = '(1I'//CNUM(LENTH)//')'
         READ (RECORD(ISTART:IEND), FMT=FORM, IOSTAT=IOS) NPI
         IF (IOS /= 0) THEN
            WRITE (6, *) 'PRSRCN: Principal quantum number ', RECORD(ISTART:&
               IEND)
            WRITE (6, *) ' could not be decoded.'
            IERR = 1
            GO TO 6
         ENDIF
         DO J = NCORE + 1, NW
            IF (NP(J)/=NPI .OR. NH(J)/=SYMI) CYCLE
            ISHELL = J
            GO TO 5
         END DO
         WRITE (6, *) 'PRSRCL: Not a peel subshell.'
         IERR = 2
         GO TO 6
    5    CONTINUE
         IOSTRT = I + 1
      ELSE IF (RECI == ')') THEN
         IOEND = I - 1
         LENTH = IOEND - IOSTRT + 1
         FORM = '(1I'//CNUM(LENTH)//')'
         READ (RECORD(IOSTRT:IOEND), FMT=FORM, IOSTAT=IOS) IOCCI
         IF (IOS /= 0) THEN
            WRITE (6, *) 'PRSRCN: Occupation number ', RECORD(IOSTRT:IOEND)
            WRITE (6, *) ' could not be decoded.'
            IERR = 3
            GO TO 6
         ENDIF
         NKJI = NKJ(ISHELL)
         IF (NKJI <= 7) THEN
            IF (IOCCI<0 .OR. IOCCI>NKJ(ISHELL)+1) THEN
               WRITE (6, *) 'PRSRCN: Occupation specified'
               WRITE (6, *) ' incorrectly for ', NP(ISHELL), NH(ISHELL)
               WRITE (6, *) ' subshell.'
               IERR = 4
               GO TO 6
            ENDIF
         ELSE
            IF (IOCCI<0 .OR. IOCCI>2) THEN
               WRITE (6, *) 'PRSRCN: Occupation specified'
               WRITE (6, *) ' incorrectly for ', NP(ISHELL), NH(ISHELL)
               WRITE (6, *) ' subshell.'
               IERR = 5
               GO TO 6
            ENDIF
         ENDIF
         IOCCS(ISHELL) = IOCCI
         ISTART = 0
      ELSE
         IF (ISTART == 0) ISTART = I
      ENDIF
!
      IF (I < 256) THEN
         I = I + 1
         GO TO 2
      ENDIF
!
      IERR = 0
!
    6 CONTINUE
      RETURN
      END SUBROUTINE PRSRCN
