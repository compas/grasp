!*******************************************************************
!                                                                  *
      SUBROUTINE RECOP00(NS,JA1,JA2,KA,IAT)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C,             ONLY: NPEEL, JJC1, JJC2, JLIST, JJQ1, JJQ2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NS, JA1, JA2, KA
      INTEGER, INTENT(OUT) :: IAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ISKR, IJ1, J, JI, JJ, NPEELGG, KK, LK, LD
!-----------------------------------------------
      IAT=1
      IF(NPEEL == 1 .AND. NS == -1)RETURN
      IAT=0
      IF(NS == -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IF(NPEELGG > 1) THEN
        LK=JJC1(NPEELGG-1)-1
        LD=JJC2(NPEELGG-1)-1
      ELSE IF(NPEELGG == 1) THEN
        LK=JJC1(NPEELGG)-1
        LD=JJC2(NPEELGG)-1
      ELSE
        PRINT*, "ERROR in RECOP00"
        STOP
      END IF
      IF(ITTK(LK,LD,2*KA) == 0)RETURN
      IAT=1
      IF(NPEEL == 1)RETURN
      DO I=1,NPEEL
        IF(JA1 /= I) THEN
          IF(JA2 /= I) THEN
            IJ1=JLIST(I)
            IF(JJQ1(1,IJ1) /= JJQ2(1,IJ1))IAT=0
            IF(JJQ1(2,IJ1) /= JJQ2(2,IJ1))IAT=0
            IF(JJQ1(3,IJ1) /= JJQ2(3,IJ1))IAT=0
          ENDIF
        ENDIF
      END DO
      IF(IAT == 0)RETURN
      IF(NPEELGG <= 2)RETURN
      IF(JA1 <= 2)RETURN
      DO J=3,JA1
        JJ=J-2
        IF(JJC1(JJ) /= JJC2(JJ))IAT=0
        IF(IAT == 0)RETURN
      END DO
      ISKR=NPEELGG-JA2
      IF(ISKR > 0) THEN
        DO JI=1,ISKR
          KK=JA2-2+JI
          LK=JJC1(KK)-1
          LD=JJC2(KK)-1
          IF(ITTK(LK,LD,2*KA) == 0)IAT=0
          IF(IAT == 0)RETURN
        END DO
      ENDIF
      RETURN
      END SUBROUTINE RECOP00
