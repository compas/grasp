!*******************************************************************
!                                                                  *
      SUBROUTINE DIAGA1(JA1,KA,IRE,IAT,RECC)
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
      INTEGER, INTENT(IN)       :: JA1, KA, IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ1,IA1,IB1,K1,J1,IT1,IT1S,IFAZ,LL1
      REAL(DOUBLE) :: A1
!-----------------------------------------------
      RECC = ZERO
      IJ1=JLIST(JA1)
      IA1=JJQ1(3,IJ1)-1
      IB1=JJQ2(3,IJ1)-1
      IF(JA1 <= 2) THEN
        LL1=JLIST(1)
        IF(JA1 == 1)LL1=JLIST(2)
        J1=JJQ1(3,LL1)-1
        IT1=JJC1(1)-1
        IT1S=JJC2(1)-1
      ELSE
        K1=JA1-2
        J1=JJC1(K1)-1
        K1=JA1-1
        IT1=JJC1(K1)-1
        IT1S=JJC2(K1)-1
      END IF
      IF(IRE /= 0) THEN
        CALL SIXJ(KA,IB1,IA1,J1,IT1,IT1S,0,A1)
        A1=A1*DSQRT(DBLE((IA1+1)*(IT1S+1)))
        IFAZ=J1+IT1+IB1+KA
        IF((IFAZ/4)*4 /= IFAZ)A1=-A1
        RECC=A1
        IAT=1
        IF(JA1 /= 1)RETURN
        IFAZ=IA1+IB1+2*J1-IT1-IT1S
        IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
      ELSE
        IF(IXJTIK(KA,IB1,IA1,J1,IT1,IT1S) == 0)RETURN
        IAT=1
      END IF
      RETURN
      END SUBROUTINE DIAGA1
