!***********************************************************************
!                                                                      *
      SUBROUTINE PRSRSL(NFILE, ID)
!                                                                      *
!   READs and parses a list of subshell labels on unit NFILE to load   *
!   COMMON  blocks  /ORB4/,  /ORB5/, and /ORB10/; the value of NW in   *
!   COMMON/ORB2/ is incremnted by 1 for each new subshell label; the   *
!   labels are delimited either by blanks or commas.                   *
!                                                                      *
!   Call(s) to: [LIB92] CONVRT.                                        *
!                                                                      *
!   Written by Farid A Parpia               Last revised 14 Sep 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:16   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE IOUNIT_C
      USE ORB_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NFILE
      INTEGER , INTENT(IN) :: ID
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISTART, I, IEND, IOS, II, IIB2, LTEGER
      LOGICAL :: NEWREC
      CHARACTER(LEN=1)   :: RECI
      CHARACTER(LEN=2)   :: CTEGER, SYM
      CHARACTER(LEN=6)   :: FORM
!      CHARACTER(LEN=500) :: RECORD, RECORD2
      CHARACTER(LEN=1000) :: RECORD, RECORD2
      CHARACTER(LEN=2), DIMENSION(19) :: SYMLST
!
      DATA SYMLST/ 's ', 'p-', 'p ', 'd-', 'd ', 'f-', 'f ', 'g-', 'g ', 'h-', &
         'h ', 'i-', 'i ', 'k-', 'k ', 'l-', 'l ', 'm-', 'm'/
!-----------------------------------------------
      IOS = 0
!
!   Read the records
!
      NEWREC = .FALSE.
      READ (NFILE, '(A)') RECORD
      IF (ID == 2) THEN
         READ (NFILE, '(A)') RECORD2
         IF (RECORD2(1:3) == 'CSF') THEN
            BACKSPACE (UNIT=NFILE)
         ELSE
            NEWREC = .TRUE.
         ENDIF
      ENDIF
!
!   Parse RECORD from left to right
!
      ISTART = 0
      I = 1
    1 CONTINUE
      RECI = RECORD(I:I)
      IF (RECI/=' ' .AND. RECI/=',') THEN
         IF (ISTART == 0) ISTART = I
      ELSE
         IF (ISTART /= 0) THEN
            IEND = I - 1
            RECI = RECORD(IEND:IEND)
            IF (RECI == '-') THEN
!              READ (RECORD(IEND-1:IEND),'(1A2)',IOSTAT = IOS) SYM
               SYM = RECORD(IEND-1:IEND)
               IF (IOS /= 0) THEN
                  WRITE (ISTDE, *) 'PRSRSL: Symmetry ', RECORD(IEND-1:IEND), &
                     ' could not be decoded.'
                  STOP
               ENDIF
               IEND = IEND - 2
            ELSE
               SYM(2:2) = ' '
!              READ (RECI,'(1A1)',IOSTAT = IOS) SYM(1:1)
               SYM(1:1) = RECI
               IF (IOS /= 0) THEN
                  WRITE (ISTDE, *) 'PRSRSL: Symmetry ', RECI, ' could not', &
                     ' be decoded.'
                  STOP
               ENDIF
               IEND = IEND - 1
            ENDIF
            DO II = 1, 19
               IF (SYM /= SYMLST(II)) CYCLE
               NW = NW + 1
               IF (NW > NNNW) THEN
                  WRITE (ISTDE, *) 'PRSRSL: Number of subshells ', &
                     'exceeds allocation: plant NW was set to', NNNW
                  STOP
               ENDIF
               NH(NW) = SYM
               IIB2 = II/2
               NKL(NW) = IIB2
               IF (MOD(II,2) == 1) THEN
                  NAK(NW) = (-IIB2) - 1
                  NKJ(NW) = II
               ELSE
                  NAK(NW) = IIB2
                  NKJ(NW) = II - 1
               ENDIF
               CALL CONVRT (IEND - ISTART + 1, CTEGER, LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND), FMT=FORM, IOSTAT=IOS) NP(NW)
               IF (IOS /= 0) THEN
                  WRITE (ISTDE, *) 'PRSRSL: Principal quantum number ', RECORD(&
                     ISTART:IEND), ' could not be decoded.'
                  STOP
               ENDIF
               GO TO 3
            END DO
            WRITE (ISTDE, *) 'PRSRSL: Symmetry ', SYM, ' could not be decoded.'
            STOP
    3       CONTINUE
            ISTART = 0
         ENDIF
      ENDIF
!
!      IF (I < 500) THEN
      IF (I < 1000) THEN
         I = I + 1
         GO TO 1
      ENDIF
      IF (NEWREC) THEN
         RECORD = RECORD2
         NEWREC = .FALSE.
         I = 1
         ISTART = 0
         GO TO 1
      ENDIF
!
      RETURN
      END SUBROUTINE PRSRSL
