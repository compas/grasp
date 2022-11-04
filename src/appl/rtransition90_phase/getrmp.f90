!***********************************************************************
!                                                                      *
      SUBROUTINE GETRMP
!                                                                      *
!   Interactively  determines the  list of radiation multipolarities   *
!   and parities. This is loadad into COMMON/OSC6/.                    *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:16:10   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE OFFD_C
      USE osc_C, ONLY: NKP, KP
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NDKP, ISTART, I, IEND, LENTH, IOS, MULT
      REAL, DIMENSION(3) :: CNUM
      LOGICAL :: LELEC, LMAGN, YES
      CHARACTER :: RECORD*256, RECI
!-----------------------------------------------
!
!
!   Initial allocation for PNTRKP
!
      NDKP = 1
      CALL ALLOC (KP, NDKP, 'KP', 'GETRMP' )
!
!   Entry message
!
    1 CONTINUE
      WRITE (6, *) 'Enter the list of transition specifications'
      WRITE (6, *) ' e.g.,  E1,M2  or  E1 M2  or  E1;M2 :'
!
!   Initialise NKP
!
    2 CONTINUE
      READ (*, '(A)') RECORD
      NKP = 0
!
!   Parse RECORD from left to right
!
      ISTART = 0
      I = 1
    3 CONTINUE
      RECI = RECORD(I:I)
      IF (RECI/=' ' .AND. RECI/=',' .AND. RECI/=';') THEN
         IF (ISTART == 0) ISTART = I
      ELSE
         IF (ISTART /= 0) THEN
            IEND = I - 1
            RECI = RECORD(ISTART:ISTART)
            IF (RECI == 'E') THEN
               LELEC = .TRUE.
               LMAGN = .FALSE.
            ELSE IF (RECI == 'M') THEN
               LELEC = .FALSE.
               LMAGN = .TRUE.
            ELSE
               WRITE (6, *) 'GETRMP: Transitions must be of type'
               WRITE (6, *) ' E or type M; reenter ...'
               GO TO 2
            ENDIF
            LENTH = IEND - ISTART
            IF (LENTH /= 1) THEN
               WRITE (6, *) 'GETRMP: Transition multipolarities'
               WRITE (6, *) ' must be integers between 1 and 9;'
               WRITE (6, *) ' reenter ...'
               GO TO 2
            ENDIF
            RECI = RECORD(IEND:IEND)
            READ (RECI, '(1I1)', IOSTAT=IOS) MULT
            IF (IOS /= 0) THEN
               WRITE (6, *) 'GETRMP: Unable to decode multipolarity'
               WRITE (6, *) ' '//RECI//'; reenter ...'
               GO TO 2
            ENDIF
            NKP = NKP + 1
            IF (NKP > NDKP) THEN
               CALL RALLOC (KP, NKP, 'KP', 'GETRMP' )
               NDKP = NKP
            ENDIF
            IF (LELEC) THEN
               KP(NKP) = MULT*(-1)**MULT
            ELSE IF (LMAGN) THEN
               KP(NKP) = MULT*(-1)**(MULT + 1)
            ENDIF
            ISTART = 0
         ENDIF
      ENDIF
!
      IF (I < 256) THEN
         I = I + 1
         GO TO 3
      ENDIF
!
      IF (NKP == 0) GO TO 1
!
!   Trim array to the exact size
!
      IF (NDKP /= NKP) CALL RALLOC (KP, NKP, 'KP', 'GETRMP')
!
!   If M1 or E2 inquire if the transitions are between levels
!   with different J quantum numbers.
!
      DO I = 1, NKP
         IF (KP(I) == 1) THEN
            WRITE (*, *) 'M1 transitions only between levels with different J?'
            YES = GETYN()
            IF (YES) THEN
               NOFFD1 = 1
            ELSE
               NOFFD1 = 0
            ENDIF
         ENDIF
         IF (KP(I) /= 2) CYCLE
         WRITE (*, *) 'E2 transitions only between levels with different J?'
         YES = GETYN()
         IF (YES) THEN
            NOFFD2 = 1
         ELSE
            NOFFD2 = 0
         ENDIF
      END DO
!
      RETURN
      END SUBROUTINE GETRMP
