!*******************************************************************
!                                                                  *
      SUBROUTINE RECO2(JA1,JA2,KA,IRE,IAT,RECC)
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
      USE m_C,             ONLY: JLIST, JJQ1, JJQ2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE diaga1_I
      USE diaga2_I
      USE diaga3_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: JA1,JA2,KA,IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IA1,IA2,IB1,IB2,IJ1,IJ2,ISKR
      REAL(DOUBLE) :: S, SS, RE
!-----------------------------------------------
      IAT=0
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      S=DBLE(JJQ1(3,IJ1))
      SS=DBLE(JJQ1(3,IJ2))
      SS=S*SS
      RECC=ONE/DSQRT(SS)
      IF(IRE == 0) THEN
        IAT=0
      ELSE IF(KA /= 0) THEN
        IAT=0
      ELSE
        IAT=1
        RETURN
      END IF
      IA1=JJQ1(3,IJ1)-1
      IA2=JJQ1(3,IJ2)-1
      IB1=JJQ2(3,IJ1)-1
      IB2=JJQ2(3,IJ2)-1
!
      CALL DIAGA2(JA1,JA2,KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC*DSQRT(DBLE(IA2+1))/DSQRT(DBLE((KA+1)*(IB2+1)))
      IF(JA1 == 1.AND.JA2 == 2)RETURN
!
      IAT=0
      CALL DIAGA1(JA1,KA,IRE,IAT,RE)
      IF(IAT == 0)RETURN
      RECC=RE*RECC
      ISKR=JA2-JA1
      IF(JA1 == 1)ISKR=JA2-1-JA1
      IF(ISKR <= 1)RETURN
!
      IAT=0
      CALL DIAGA3(JA1,JA2,KA,IRE,IAT,RE)
      RECC=RE*RECC
      RETURN
      END SUBROUTINE RECO2
