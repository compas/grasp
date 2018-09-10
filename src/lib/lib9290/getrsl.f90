!***********************************************************************
!                                                                      *
      SUBROUTINE GETRSL(INDX, NSUBS) 
!                                                                      *
!   READs and parses a list of relativistic subshell labels delimit-   *
!   ed either by blanks or by commas. An asterisk may be used as the   *
!   `wildcard' character.                                              *
!
!   Output:
!     NSUBS - the # of orbitals parsed
!     indx(1:NSUBS) - indeces of these orbitals
!                                                                      *
!   Call(s) to: [LIB92]: LDIGIT.                                       *
!                                                                      *
!   Written by Farid A. Parpia             Last revised: 18 Dec 1992   *
!   Modified by Xinghong He                Last revised: 09 Jul 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:14   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE ORB_C 
      USE IOUNIT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ldigit_I 
      USE convrt_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: NSUBS 
      INTEGER, DIMENSION(*), INTENT(INOUT) :: INDX
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISTART, I, IEND, NFORM, J, LENTH, IOS, NPQ, K 
      LOGICAL :: FOUND, NLANY, NLOK, NPANY, NPOK, NSANY, NSOK 
      CHARACTER :: NLQ, NSQ, RECI, CNUM*2, FORM*5, RECORD*500 
!-----------------------------------------------
!
!
 
      GO TO 2 
 
    1 CONTINUE 
      WRITE (ISTDE, *) ' redo ...' 
!
!   Read a record
!
    2 CONTINUE 
      NSUBS = 0 
      READ (*, '(A)') RECORD 

      WRITE(734,'(a)') trim(record)
!
!   Parse RECORD from left to right
!
      ISTART = 0 
      I = 1 
    3 CONTINUE 
      RECI = RECORD(I:I) 
      IF (RECI/=' ' .AND. RECI/=',') THEN 
         IF (ISTART == 0) ISTART = I 
      ELSE 
         IF (ISTART /= 0) THEN 
            IEND = I - 1 
!
!   Parse the substring from left to right
!
!   (1) Determine the principal quantum number
!
            IF (RECORD(ISTART:ISTART) == '*') THEN 
               NPANY = .TRUE. 
               ISTART = MIN(ISTART + 1,IEND) 
            ELSE 
               NPANY = .FALSE. 
               NFORM = 0 
               DO J = ISTART, IEND 
                  IF (.NOT.LDIGIT(RECORD(J:J))) CYCLE  
                  NFORM = NFORM + 1 
               END DO 
               IF (NFORM == 0) THEN 
                  WRITE (ISTDE, *) 'GETRSL: Unable to interpret ', &
                     'the principal quantum number;' 
                  GO TO 1 
               ENDIF 
               CALL CONVRT (NFORM, CNUM, LENTH) 
               FORM = '(1I'//CNUM(1:LENTH)//')' 
               READ (RECORD(ISTART:ISTART+NFORM-1), FORM, IOSTAT=IOS) NPQ 
               IF (IOS /= 0) THEN 
                  WRITE (ISTDE, *) 'GETRSL: Unable to interpret ', 'string ', &
                     RECORD(ISTART:IEND-2), ' as a principal quantum number' 
                  GO TO 1 
               ENDIF 
               ISTART = ISTART + NFORM 
            ENDIF 
!
!   (2) Determine the orbital angular momentum quantum number
!
            NLQ = RECORD(ISTART:ISTART) 
            IF (NLQ == '*') THEN 
               NLANY = .TRUE. 
            ELSE 
               NLANY = .FALSE. 
            ENDIF 
!
!   (3) Determine the spin-orbit component
!
            IF (IEND > ISTART) THEN 
               NSQ = RECORD(IEND:IEND) 
               IF (NSQ == '*') THEN 
                  NSANY = .TRUE. 
               ELSE IF (NSQ == '-') THEN 
                  NSANY = .FALSE. 
               ELSE 
                  WRITE (ISTDE, *) 'GETRSL: Unable to interpret ', 'string ', &
                     NSQ, ' as a spin-orbit component indicator' 
                  GO TO 1 
               ENDIF 
            ELSE 
               IF (NLANY) THEN 
                  NSANY = .TRUE. 
               ELSE 
                  NSANY = .FALSE. 
                  NSQ = ' ' 
               ENDIF 
            ENDIF 
!
!
            FOUND = .FALSE. 
            DO J = 1, NW 
               NPOK = NPANY .OR. NP(J)==NPQ 
               NLOK = NLANY .OR. NLQ==NH(J)(1:1) 
               NSOK = NSANY .OR. NSQ==NH(J)(2:2) 
               IF (.NOT.(NPOK .AND. NLOK .AND. NSOK)) CYCLE  
               DO K = 1, NSUBS 
                  IF (INDX(K) /= J) CYCLE  
                  WRITE (ISTDE, *) 'GETRSL: ', 'Repeated subshell in list;' 
                  GO TO 1 
               END DO 
               FOUND = .TRUE. 
               NSUBS = NSUBS + 1 
               INDX(NSUBS) = J 
            END DO 
!
            IF (.NOT.FOUND) THEN 
               WRITE (ISTDE, *) 'GETRSL: Subshell not occupied as ', &
                  ' according to CSL  File;' 
               GO TO 1 
            ENDIF 
!
            ISTART = 0 
         ENDIF 
      ENDIF 
!
      IF (I < 500) THEN 
         I = I + 1 
         GO TO 3 
      ENDIF 
!
      RETURN  
      END SUBROUTINE GETRSL 
