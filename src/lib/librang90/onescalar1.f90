!*******************************************************************
!                                                                  *
      SUBROUTINE ONESCALAR1(NS,JJA,JJB,JA,JB,COEFF)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
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
      USE CONS_C,          ONLY: ZERO, HALF, ONE, EPS
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoonescalar_I
      USE perko2_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NS,JJA,JJB,JA,JB
      REAL(DOUBLE), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IAT
      REAL(DOUBLE) :: WJ,QM1,QM2,RECOUPL
!-----------------------------------------------
!
!     THE CASE 1111   + + - -
!
      COEFF=ZERO
      IF(JA == JB) THEN
        IF(JJA /= JJB) THEN
          CALL RECOONESCALAR(NS,JA,JA,JA,JA,0,IAT)
          IF(IAT == 0)RETURN
        END IF
        CALL PERKO2(JA,JA,JA,JA,1)
        QM1=HALF
        QM2=-HALF
        CALL WJ1(IK1,BK1,ID1,BD1,0,QM1,QM2,WJ)
        IF(DABS(WJ) > EPS) THEN
           RECOUPL=ONE/DSQRT(DBLE(IK1(6)+1))
           COEFF=WJ*RECOUPL*DSQRT(DBLE(ID1(3)+1))
           COEFF=-COEFF
        END IF
      END IF
      RETURN
      END SUBROUTINE ONESCALAR1
