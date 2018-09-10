!*******************************************************************
!                                                                  *
      SUBROUTINE RECO4(JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 09  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
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
      USE m_C,             ONLY: JLIST, JJQ1, JJQ2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE diaga1_I
      USE diaga2_I
      USE diaga3_I
      USE diaga4_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,IRE
      INTEGER, INTENT(OUT) :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IA4, IB4, IJ1, IJ2, IJ3, IJ4, ISKR
      REAL(DOUBLE) :: S, S1, S2, S3, S4, RE
!-----------------------------------------------
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IJ3=JLIST(JA3)
      IJ4=JLIST(JA4)
      S1=JJQ1(3,IJ1)
      S2=JJQ1(3,IJ2)
      S3=JJQ1(3,IJ3)
      S4=JJQ1(3,IJ4)
      S=S1*S2*S3*S4
      RECC=ONE/DSQRT(S)
      IA4=JJQ1(3,IJ4)-1
      IB4=JJQ2(3,IJ4)-1
      RECC=RECC*DSQRT(DBLE(IA4+1))/DSQRT(DBLE((K4+1)*(IB4+1)))
!
      ISKR=JA3-JA2
      IF(ISKR > 1) THEN
        IAT=0
        CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
        IF(IAT == 0)RETURN
        RECC=RE*RECC
      END IF
!
      ISKR=JA4-JA3
      IF(ISKR > 1) THEN
        IAT=0
        CALL DIAGA3(JA3,JA4,K4,IRE,IAT,RE)
        IF(IAT == 0)RETURN
        RECC=RE*RECC
      END IF
!
      IAT=0
      CALL DIAGA2(JA1,JA4,K4,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
!
      IAT=0
      CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
!
      IAT=0
      CALL DIAGA4(JA2,JA3,KA,K3,K4,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
      IF(JA1 == 1.AND.JA2 == 2)RETURN
!
      IAT=0
      CALL DIAGA1(JA1,K1,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
!
      ISKR=JA2-JA1
      IF(JA1 == 1)ISKR=JA2-1-JA1
      IF(ISKR <= 1)RETURN
      IAT=0
      CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
      RECC=RE*RECC
      RETURN
      END SUBROUTINE RECO4
