!***********************************************************************
!                                                                      *
      SUBROUTINE SETQNA(JA, JB) 
!                                                                      *
!   This generates the  arrays  defining  the quantum numbers of the   *
!   states involved in the  matrix  element  linking  configurations   *
!   labelled by JA, JB.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ICHOP, IQ, JCUP, JQS.                         *
!                                                                      *
!                                           Last update: 30 Oct 1987   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:39   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE M_C
      USE ORB_C,           ONLY: NCF, NW 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I 
      USE jqs_I 
      USE ichop_I 
      USE jcup_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JA 
      INTEGER  :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, JCNT, JCNTOP, JW1, JW2, JW 
!-----------------------------------------------
!
!
!   List parameters defining all shells in both configurations, whether
!   participating or not
!
      DO J = 1, NW 
         NQ1(J) = IQ(J,JA) 
         NQ2(J) = IQ(J,JB) 
         DO K = 1, 3 
            JJQ1(K,J) = JQS(K,J,JA) 
            JJQ2(K,J) = JQS(K,J,JB) 
         END DO 
      END DO 
!
!   Define coupling schemes: set JLIST array to define those shells
!   which are open in either configuration, and KLIST array to locate
!   the rest. Exclude shells which are empty in both configurations
!
      NPEEL = 0 
      NCORE = 0 
      DO J = 1, NW 
         IF (ICHOP(J,JA)==(-1) .AND. ICHOP(J,JB)==(-1)) CYCLE  
         IF (ICHOP(J,JA)==1 .AND. ICHOP(J,JB)==1) THEN 
            NCORE = NCORE + 1 
            KLIST(NCORE) = J 
         ELSE 
            NPEEL = NPEEL + 1 
            JLIST(NPEEL) = J 
         ENDIF 
      END DO 
!
!   Return if not more than one shell is open
!
      IF (NPEEL <= 1) RETURN  
!
!   Set arrays of coupling angular momenta interpolating closed
!   shells where necessary. Left hand side first ...
!
      JCNT = 1 
      JCNTOP = 0 
      JW1 = JLIST(1) 
      JW2 = JLIST(2) 
      IF (ICHOP(JW1,JA) /= 0) THEN 
         JJC1(1) = JQS(3,JW2,JA) 
         IF (ICHOP(JW2,JA) == 0) JCNTOP = 1 
      ELSE 
         JCNTOP = 1 
         IF (ICHOP(JW2,JA) == 0) THEN 
            JJC1(1) = JCUP(JCNT,JA) 
            JCNT = JCNT + 1 
         ELSE 
            JJC1(1) = JQS(3,JW1,JA) 
         ENDIF 
      ENDIF 
!
      DO J = 3, NPEEL 
         JW = JLIST(J) 
         IF (ICHOP(JW,JA) /= 0) THEN 
            JJC1(J-1) = JJC1(J-2) 
         ELSE 
            IF (JCNTOP /= 0) THEN 
               JJC1(J-1) = JCUP(JCNT,JA) 
               JCNT = JCNT + 1 
            ELSE 
               JJC1(J-1) = JQS(3,JW,JA) 
            ENDIF 
            JCNTOP = JCNTOP + 1 
         ENDIF 
      END DO 
!
!   ... and repeat for right hand side
!
      JCNT = 1 
      JCNTOP = 0 
      JW1 = JLIST(1) 
      JW2 = JLIST(2) 
      IF (ICHOP(JW1,JB) /= 0) THEN 
         JJC2(1) = JQS(3,JW2,JB) 
         IF (ICHOP(JW2,JB) == 0) JCNTOP = 1 
      ELSE 
         JCNTOP = 1 
         IF (ICHOP(JW2,JB) == 0) THEN 
            JJC2(1) = JCUP(JCNT,JB) 
            JCNT = JCNT + 1 
         ELSE 
            JJC2(1) = JQS(3,JW1,JB) 
         ENDIF 
      ENDIF 
!
      DO J = 3, NPEEL 
         JW = JLIST(J) 
         IF (ICHOP(JW,JB) /= 0) THEN 
            JJC2(J-1) = JJC2(J-2) 
         ELSE 
            IF (JCNTOP /= 0) THEN 
               JJC2(J-1) = JCUP(JCNT,JB) 
               JCNT = JCNT + 1 
            ELSE 
               JJC2(J-1) = JQS(3,JW,JB) 
            ENDIF 
            JCNTOP = JCNTOP + 1 
         ENDIF 
      END DO 
!
      RETURN  
      END SUBROUTINE SETQNA 
