!*******************************************************************
!                                                                  *
      SUBROUTINE RECOP2(NS,JA1,JA2,K1,K2,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ONE
      USE m_C,             ONLY: JLIST, JJQ1, JJQ2, NPEEL
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE diaga1_I
      USE diaga3_I
      USE diaga4_I
      USE diaga5_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NS, JA1, JA2, K1, K2, KA, IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ, ISKR, NPEELGG
      REAL(DOUBLE) :: S, SS, RE
!-----------------------------------------------
      IAT=1
      IJ=JLIST(JA1)
      S=DBLE(JJQ1(3,IJ))
      IJ=JLIST(JA2)
      SS=DBLE(JJQ1(3,IJ))
      SS=S*SS
      RECC=ONE/DSQRT(SS)
      IF(NPEEL == 1 .AND. NS == -1)RETURN
      IF(NS == -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IAT=0
      ISKR=NPEELGG-JA2
      IF(ISKR > 1) THEN
        CALL DIAGA3(JA2,NPEELGG,2*KA,IRE,IAT,RE)
        IF(IAT == 0)RETURN
        RECC=RE*RECC
        IAT=0
      END IF
      IF(JA2 /= NPEELGG) THEN
        CALL DIAGA5(NPEELGG,JA2,2*KA,IRE,IAT,RE)
        IF(IAT == 0)RETURN
        RECC=RE*RECC
        IAT=0
      ENDIF
      CALL DIAGA4(JA1,JA2,K1,K2,2*KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
      IF(JA1 == 1.AND.JA2 == 2)RETURN
      IAT=0
      CALL DIAGA1(JA1,K1,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
      ISKR=JA2-JA1
      IF(JA1 == 1)ISKR=JA2-1-JA1
      IF(ISKR <= 1)RETURN
      IAT=0
      CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
      RECC=RE*RECC
      RETURN
      END SUBROUTINE RECOP2
