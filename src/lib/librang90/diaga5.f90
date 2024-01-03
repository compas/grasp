!*******************************************************************
!                                                                  *
      SUBROUTINE DIAGA5(NPEELGG,JA1,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 01  --------------  *
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
      USE m_C,             ONLY: JLIST, JJC1, JJC2, JJQ1, JJQ2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NPEELGG,JA1,KA,IRE
!      INTEGER, INTENT(OUT)      :: IAT
      INTEGER, INTENT(INOUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: ITI1,ITI1S,IJ1,ITI,ITIS,JI
      REAL(DOUBLE) :: A3
!-----------------------------------------------
      RECC = ZERO
      ITI1=JJC1(NPEELGG-1)-1
      ITI1S=JJC2(NPEELGG-1)-1
      IJ1=JLIST(NPEELGG)
      IF(JA1 == NPEELGG) THEN
        ITI=JJQ1(3,IJ1)-1
        ITIS=JJQ2(3,IJ1)-1
        JI=JJC1(NPEELGG-2)-1
      ELSE
        JI=JJQ1(3,IJ1)-1
        ITI=JJC1(NPEELGG-2)-1
        ITIS=JJC2(NPEELGG-2)-1
      END IF
      IF(IRE == 0) THEN
        IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) /= 0) IAT=1
      ELSE
        CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
        RECC=A3*DSQRT(DBLE((ITI+1)*(ITI1S+1)))
        IF(MOD(KA+JI+ITIS+ITI1,4) /= 0)RECC=-RECC
        IAT=1
        IF(JA1 == NPEELGG)RETURN
        IF(MOD(ITI+ITIS-ITI1S-ITI1+2*JI,4) /= 0)RECC=-RECC
      END IF
      RETURN
      END SUBROUTINE DIAGA5
