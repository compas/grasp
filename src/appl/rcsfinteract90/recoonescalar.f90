!*******************************************************************
!                                                                  *
      SUBROUTINE RECOONESCALAR(NS,JA1,JA2,JA3,JA4,KA,IAT)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NS, JA1, JA2, JA3, JA4, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IA1, IA2, IJ, IJ1, IJ2, J, NPEELGG
!-----------------------------------------------
      IAT=1
      IF(NPEEL.EQ.1 .AND. NS.EQ.-1)RETURN
      IF(NS .EQ. -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1.EQ.1.AND.JA2.EQ.2) GO TO 1
      IF(KA.NE.0)GO TO 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I.LT.NPEELGG-1) THEN
         IF(JJC1(I).NE.JJC2(I))IAT=0
        END IF
        IF(KA.NE.0) THEN
          IF(I.EQ.JA1) CYCLE
          IF(I.EQ.JA2) CYCLE
        END IF
        DO J=1,3
          IF(JJQ1(J,IJ).NE.JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
!
!  OTHER CASES
!
    5 CONTINUE
      DO I=1,NPEELGG
        IJ=JLIST(I)
        IF(I.LT.NPEELGG-1) THEN
          IA1=JA1-1
          IA2=JA2-1
          IF(JA1.EQ.1)IA1=JA1
          IF(I.GE.IA1.AND.I.LT.IA2)GO TO 7
          IF(JJC1(I).NE.JJC2(I))IAT=0
        END IF
    7   IF(I.EQ.JA1) CYCLE
        IF(I.EQ.JA2) CYCLE
        IF((KA.EQ.2).AND.(I.EQ.JA3)) CYCLE
        IF((KA.EQ.3).AND.(I.EQ.JA3)) CYCLE
        IF((KA.EQ.3).AND.(I.EQ.JA4)) CYCLE
        DO J=1,3
          IF(JJQ1(J,IJ).NE.JJQ2(J,IJ))IAT=0
        END DO
      END DO
      RETURN
      END SUBROUTINE RECOONESCALAR
