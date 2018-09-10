!*******************************************************************
!                                                                  *
      SUBROUTINE ONEPARTICLEJJ1(NS,KA,JJA,JJB,JA,JB,COEFF)
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
      USE CONS_C,          ONLY: ZERO, HALF, EPS
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recop00_I
      USE recop1_I
      USE recop2_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NS,KA,JJA,JJB,JA,JB
      REAL(DOUBLE), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IAT
      REAL(DOUBLE) :: REC,WJ,QM1,QM2
!-----------------------------------------------
!
!     THE CASE 11   + -
!
      COEFF=ZERO
      IF(JA == JB) THEN
        IF(JJA /= JJB) THEN
          CALL RECOP00(NS,JA,JA,KA,IAT)
          IF(IAT == 0)RETURN
        END IF
        CALL RECOP1(NS,JA,KA,0,IAT,REC)
        IF(IAT == 0)RETURN
        CALL PERKO2(JA,JA,JA,JA,1)
        QM1=HALF
        QM2=-HALF
        CALL WJ1(IK1,BK1,ID1,BD1,KA,QM1,QM2,WJ)
        IF(DABS(WJ) > EPS) THEN
           CALL RECOP1(NS,JA,KA,1,IAT,REC)
           COEFF=WJ*REC*DSQRT(DBLE(ID1(3)+1))/DSQRT(DBLE(2*KA+1))
           COEFF=-COEFF
        END IF
      END IF
      RETURN
      END SUBROUTINE ONEPARTICLEJJ1
