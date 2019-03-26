!***********************************************************************
!                                                                      *
      SUBROUTINE PARSJL(MODE, NCORE, RECORD, LOC, JX, NJX, IERR)
!                                                                      *
!   READs and  parses a string that specifies angular momentum quan-   *
!   tum numbers.                                                       *
!                                                                      *
!   Call(s) to: [LIB92] CONVRT.                                        *
!                                                                      *
!   Written by Farid A Parpia              Last revised: 21 Dec 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:09   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      USE ORB_C, ONLY: NCF, NW, IQA
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)           :: MODE, NCORE, LOC
      INTEGER, INTENT(OUT)          :: NJX, IERR
      CHARACTER(LEN=256), INTENT(IN) :: RECORD
      INTEGER, DIMENSION(*), INTENT(OUT) :: JX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NJXMAX, ISTART, I, ILOOP, IFRAC, IEND, LTEGER, J
      CHARACTER(LEN=1) :: RECI
      CHARACTER(LEN=2) :: CTEGER
      CHARACTER(LEN=6) :: FORM
!-----------------------------------------------
!
! This subroutine is here to replace the existing one which has been
! renamed as PARSJL_OLD. The purpose is to remove illegal GOTO into an
! IF block. The strategy is simple: copy the kernl to the bottom
! and add a do-loop over I.
! To restore (re-use) this subroutine, give the existing PARSJL a new
! name and then change the name of this subroutine back to PARSJL, and
! compile.
! XHH 1997.01.28
!
!
!   There cannot be more than JXMAX angular momenta specified; if
!   MODE is 1, the subshell quantum numbers are being read; if MODE
!   is 2, the intermediate and final angular momenta are being read
!
      IF (MODE == 1) THEN
         NJXMAX = 2*(NW - NCORE)
      ELSE
         NJXMAX = NW - NCORE
      ENDIF
!
!   Initialise NJX
!
      NJX = 0
!
!   Parse RECORD from left to right
!
      ISTART = 0
      I = 1

! The original algorithm goes through the whole subroutine at least
! once, whatever the value of LOC. Thus we define another integer
! iloop to achieve this.
! XHH 1997.01.28

      ILOOP = MAX(1,LOC)
      DO I = 1, ILOOP
         RECI = RECORD(I:I)
         IF (RECI/=' ' .AND. RECI/=',' .AND. RECI/=';') THEN
            IF (ISTART == 0) THEN
               ISTART = I
               IFRAC = 0
            ELSE
               IF (RECI == '/') THEN
                 IFRAC = I
               ENDIF
            ENDIF
         ELSE
            IF (ISTART /= 0) THEN
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               NJX = NJX + 1
               IF (NJX > NJXMAX) THEN
                  WRITE (6, *) 'PARSJL: Too many angular momentum', &
                     ' quantum numbers specified;'
                  IERR = 1
                  GO TO 3
               ENDIF
               IEND = I - 1
               IF (IFRAC == 0) THEN
                  CALL CONVRT (IEND - ISTART + 1, CTEGER, LTEGER)
                  FORM = '(1I'//CTEGER(1:LTEGER)//')'
                  READ (RECORD(ISTART:IEND), FMT=FORM) J
                  IF (J < 0) THEN
                     WRITE (6, *) 'PARSJL: Negative angular momentum', &
                        ' quantum number found;'
                     IERR = 2
                     GO TO 3
                  ENDIF
                  JX(NJX) = 2*J
               ELSE
                  CALL CONVRT (IEND - IFRAC, CTEGER, LTEGER)
                  FORM = '(1I'//CTEGER(1:LTEGER)//')'
                  READ (RECORD(IFRAC+1:IEND), FMT=FORM) J
                  IF (J /= 2) THEN
                     WRITE (6, *) 'PARSJL: The denominator of a', &
                        ' fractional quantum number must be 2;'
                     IERR = 3
                     GO TO 3
                  ENDIF
                  CALL CONVRT (IFRAC - ISTART, CTEGER, LTEGER)
                  FORM = '(1I'//CTEGER(1:LTEGER)//')'
                  READ (RECORD(ISTART:IFRAC-1), FMT=FORM) J
                  IF (J < 0) THEN
                     WRITE (6, *) 'PARSJL: Negative angular momentum', &
                        ' quantum number found;'
                     IERR = 4
                     GO TO 3
                  ENDIF
                  JX(NJX) = J
               ENDIF
               ISTART = 0
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ENDIF
         ENDIF
      END DO
!
      IF (LOC >= 1) THEN

! The following was accessed one extra time only when 1 <= I = LOC+1 .
! After the do-loop above, the value of I would be either LOC+1 or
! 2, depending on if LOC is greater than 1 or not, respectively


!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The following is exactly the same as those above
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         NJX = NJX + 1
         IF (NJX > NJXMAX) THEN
            WRITE (6, *) 'PARSJL: Too many angular momentum', &
               ' quantum numbers specified;'
            IERR = 1
            GO TO 3
         ENDIF
         IEND = I - 1
         IF (IFRAC == 0) THEN
            CALL CONVRT (IEND - ISTART + 1, CTEGER, LTEGER)
            FORM = '(1I'//CTEGER(1:LTEGER)//')'
            READ (RECORD(ISTART:IEND), FMT=FORM) J
            IF (J < 0) THEN
               WRITE (6, *) 'PARSJL: Negative angular momentum', &
                  ' quantum number found;'
               IERR = 2
               GO TO 3
            ENDIF
            JX(NJX) = 2*J
         ELSE
            CALL CONVRT (IEND - IFRAC, CTEGER, LTEGER)
            FORM = '(1I'//CTEGER(1:LTEGER)//')'
            READ (RECORD(IFRAC+1:IEND), FMT=FORM) J
            IF (J /= 2) THEN
               WRITE (6, *) 'PARSJL: The denominator of a', &
                  ' fractional quantum number must be 2;'
               IERR = 3
               GO TO 3
            ENDIF
            CALL CONVRT (IFRAC - ISTART, CTEGER, LTEGER)
            FORM = '(1I'//CTEGER(1:LTEGER)//')'
            READ (RECORD(ISTART:IFRAC-1), FMT=FORM) J
            IF (J < 0) THEN
               WRITE (6, *) 'PARSJL: Negative angular momentum', &
                  ' quantum number found;'
               IERR = 4
               GO TO 3
            ENDIF
            JX(NJX) = J
         ENDIF
         ISTART = 0
      ENDIF
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IERR = 0
!
    3 CONTINUE
      RETURN
      END SUBROUTINE PARSJL
