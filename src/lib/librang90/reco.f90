!*******************************************************************
!                                                                  *
      SUBROUTINE RECO(JA1,JA2,JA3,JA4,KA,IAT)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2020  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C,             ONLY: NPEEL, JLIST, JJC1, JJC2, JJQ1, JJQ2
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA1, JA2, JA3, JA4, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IJ1, IJ2, IJ, J, IA1, IA2
!-----------------------------------------------
      IAT=1
      IF(NPEEL == 1)RETURN
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1 == 1.AND.JA2 == 2) GO TO 1
      IF(KA /= 0)GO TO 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 DO I=1,NPEEL
        IJ=JLIST(I)
        IF(I < NPEEL-1) THEN
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
        IF(KA /= 0) THEN
          IF(I == JA1) CYCLE
          IF(I == JA2) CYCLE
        END IF
        DO J=1,3
          IF(JA1 == JA2 .AND. I == JA1 .AND. J == 1) CYCLE
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
!
!  OTHER CASES
!
    5 CONTINUE
      DO I=1,NPEEL
        IJ=JLIST(I)
        IF(I < NPEEL-1) THEN
          IA1=JA1-1
          IA2=JA2-1
          IF(JA1 == 1)IA1=JA1
          IF(I >= IA1.AND.I < IA2)GO TO 7
          IF(JJC1(I) /= JJC2(I))IAT=0
        END IF
    7   IF(I == JA1) CYCLE
        IF(I == JA2) CYCLE
        IF((KA == 2).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA3)) CYCLE
        IF((KA == 3).AND.(I == JA4)) CYCLE
        DO J=1,3
          IF(JJQ1(J,IJ) /= JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
      END SUBROUTINE RECO
