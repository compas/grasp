!*******************************************************************
!                                                                  *
      SUBROUTINE DIAGA2(JA1,JA2,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 02  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
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
      USE ixjtik_I
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: JA1, JA2, KA, IRE
!      INTEGER, INTENT(OUT)      :: IAT
      INTEGER, INTENT(INOUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ1,IJ2,IA1,IA2,IB1,IB2,N1,N2,J2,IT2,IT2S,IFAZ
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
      ELSE
        N1=JA2-1
        J2=JJC1(N1)-1
        N2=JA2-2
        IT2=JJC1(N2)-1
        IT2S=JJC2(N2)-1
      END IF
      IF(IRE /= 0) THEN
        CALL SIXJ(KA,IB2,IA2,J2,IT2,IT2S,0,A2)
        RECC=A2*DSQRT(DBLE((IB2+1)*(IT2+1)))
        IFAZ=J2+IT2S+IA2+KA
        IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
        IAT=1
      ELSE
        IF(IXJTIK(KA,IB2,IA2,J2,IT2,IT2S) == 0)RETURN
        IAT=1
      END IF
      RETURN
      END SUBROUTINE DIAGA2
