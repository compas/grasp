!*******************************************************************
!                                                                  *
      SUBROUTINE ONESCALARINT(JA,JB,IA1,IA2,INTERACT)
!                                                                  *
!   The main program for evaluating the reduced matrix elements of *
!   a one particle operator for configurations in jj-coupling.     *
!                                                                  *
!   Call(s) to: [LIB92]: CFP, FIXJ, GENSUM, ICHOP, IROW1, ISPAR,   *
!                        ITJPO, ITRIG, SETQNA, VIJOUT.             *
!                                                                  *
!   Written by  G. Gaigalas                   NIST, December 2015  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE m_C
      USE orb_C,   ONLY: NW, NAK
      USE dumx_C,  ONLY: JLIS, JC1S, JC2S
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itrig_I
      USE ichkq1_I
      USE ichop_I
      USE ispar_I
      USE itjpo_I
      USE setqna_I
      USE onescalar1INT_I
      USE onescalar2INT_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA,JB
      INTEGER, INTENT(OUT) :: IA1,IA2
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK,IOPAR,IJ,IDQ,JA1,JA2,NDQ,NS,ISH,I,II,I1,IM, &
                 JW,KS1,KS2,NX,NPEELM
      INTEGER, DIMENSION(2) :: IS,KS
!-----------------------------------------------
      IA1 = 0
      KK = 1
      IOPAR = 1
      INTERACT = 0
      IF(ITRIG(ITJPO (JA),ITJPO (JB),KK).EQ.0)RETURN
      IF((IOPAR.NE.0).AND.(ISPAR(JA)*ISPAR(JB)*IOPAR.NE.1))RETURN
      IF(ICHKQ1(JA,JB).EQ.0)RETURN
!
      CALL SETQNA (JA,JB)
!
!
!   Analyse peel shell interactions
      IDQ = 0
      JA1 = 0
      JA2 = 0
      IF (NPEEL .NE. 0) THEN
        DO JW = 1,NPEEL
          IJ = JLIST(JW)
          NDQ = NQ1(IJ)-NQ2(IJ)
          IF (IABS (NDQ) .GT. 1) RETURN
          IF (NDQ .GT. 0) THEN
            JA1 = JW
            IDQ = IDQ+1
          ELSEIF (NDQ .LT. 0) THEN
            JA2 = JW
            IDQ = IDQ+1
          ENDIF
        END DO
!
        IF (IDQ .GT. 2) RETURN
!
!   Evaluate the array VSHELL
!
!   Then there are two possibilities IDQ = 0 or IDQ = 2
!   if IDQ = 0, then loop over all shells by index ISH
!   if IDQ = 2, then one orbital fixed on each side
        NS = NPEEL
      ENDIF
!
      IF (IDQ .EQ. 2) GOTO 19
!
!   Loop over shells when IDQ = 0
      ISH = 0
      IF (NPEEL .EQ. 0) GOTO 9
      DO I = 1,NPEEL
         JLIS(I) = JLIST(I)
      END DO
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO I = 1,NPEELM
         JC1S(I) = JJC1(I)
         JC2S(I) = JJC2(I)
      END DO
!
!   If ISH .GT. NW, then loop is over and return
    9 ISH = ISH+1
      IF (ISH .GT. NW) RETURN
      IF (ICHOP (ISH,JA) .EQ. -1) GOTO 9
      IF (ICHOP (ISH,JA) .EQ. 0) GOTO 16
!
!   Case one --- the ISH-th shell is in the core or in the peel and
!   closed for both sides
      I = 1
      IF (NPEEL.EQ.0) GOTO 15
      DO I = 1,NPEEL
        IJ = JLIST(I)
        IF (ISH .LT. IJ) GOTO 11
      END DO
      I = NPEEL+1
      GOTO 13
   11 IM = NPEEL-I+1
      DO II = 1,IM
         JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
         IF (NPEEL.EQ.II) GOTO 13
         JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
         JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
      END DO
   13 CONTINUE
      IF (I .LT. 3) GOTO 14
      JJC1(I-1) = JJC1(I-2)
      JJC2(I-1) = JJC2(I-2)
      GOTO 15
   14 I1 = JLIST(1)
      JJC1(1) = JJQ1(3,I1)
      JJC2(1) = JJQ2(3,I1)
   15 JLIST(I) = ISH
      JA1 = I
      JA2 = I
      NS = NPEEL+1
      GOTO 19
!
!   Case two --- the ISH-th shell is in the peel and open for either
!   side
   16 NS = NPEEL
      DO  JW = 1,NPEEL
        NX = ISH-JLIST(JW)
        IF (NX.EQ.0) GOTO 18
      END DO
   18 JA1 = JW
      JA2 = JW
!
!   Main computation
!
!     JA1, JA2 are the indices of interacting shells in JLIST
!     IA1, IA2 are the indices of interacting shells in NW
   19 IA1 = JLIST(JA1)
      IA2 = JLIST(JA2)
      KS1 = 2*IABS (NAK(IA1))
      KS2 = 2*IABS (NAK(IA2))
!
!   Check triangular condition for the active shells
      IF (ITRIG (KS1,KS2,KK).EQ.1) GOTO 99
      IF (IDQ .EQ. 2) RETURN
      GOTO 100
!
!   Set tables of quantum numbers of non-interacting spectator shells
   99 CONTINUE

      IF (IDQ .EQ. 0) THEN
!GG
!GG   Cia orginalioje programoje yra klaida
!GG
        IF(JA .EQ. JB) THEN
          CALL ONESCALAR1INT(NS,JA,JB,JA1,JA2,INTERACT)
        ELSE
          INTERACT = 0
        END IF
      ELSE IF (IDQ .EQ. 2) THEN
!
!   IDQ = 2 Case
!
!       Permutation factor for IDQ = 2
        CALL ONESCALAR2INT(JA,JB,JA1,JA2,INTERACT)
        RETURN
      END IF
!
!   End of loop over parent states
!
!
!   IDQ = 0 CASE
!
!   Loop over all shells when IDQ = 0
100   CONTINUE
      IF (NPEEL .EQ. 0) GOTO 9
      DO I = 1,NPEEL
         JLIST(I) = JLIS(I)
      END DO
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO  I = 1,NPEELM
         JJC1(I)  = JC1S(I)
         JJC2(I)  = JC2S(I)
      END DO
      GOTO 9
      RETURN
      END SUBROUTINE ONESCALARINT
