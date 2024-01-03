!*******************************************************************
!                                                                  *
      SUBROUTINE DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 04  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  NINE                                     *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO
      USE m_C,             ONLY: JLIST, JJQ1, JJQ2, JJC1, JJC2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nine_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: JA1,JA2,K1,K2,KA,IRE
!      INTEGER, INTENT(OUT)      :: IAT
      INTEGER, INTENT(INOUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,J2S,IT2,IT2S
      REAL(DOUBLE) :: A2
!-----------------------------------------------
      RECC = ZERO
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IA1=JJQ1(3,IJ1)-1
      IA2=JJQ1(3,IJ2)-1
      IB1=JJQ2(3,IJ1)-1
      IB2=JJQ2(3,IJ2)-1
      IF(JA1 == 1.AND.JA2 == 2) THEN
        IT2=IA1
        IT2S=IB1
        J2=JJC1(1)-1
        J2S=JJC2(1)-1
      ELSE
        N1=JA2-1
        J2=JJC1(N1)-1
        J2S=JJC2(N1)-1
        N2=JA2-2
        IT2=JJC1(N2)-1
        IT2S=JJC2(N2)-1
      END IF
      IF(IRE /= 0) THEN
        CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,0,IAT,A2)
        RECC=A2*DSQRT(DBLE((IT2+1)*(KA+1)*(IA2+1)*(J2S+1)))
      ELSE
        CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,1,IAT,A2)
      END IF
      RETURN
      END SUBROUTINE DIAGA4
