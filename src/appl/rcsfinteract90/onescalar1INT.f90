!*******************************************************************
!                                                                  *
      SUBROUTINE ONESCALAR1INT(NS,JJA,JJB,JA,JB,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
!                                                                  *
!   Written by  G. Gaigalas                   NIST, December 2015  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
!      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoonescalar_I
!      USE perko2_I
!      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NS,JJA,JJB,JA,JB
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IAT
!      REAL(DOUBLE) :: WJ,QM1,QM2,RECOUPL
!-----------------------------------------------
!
!     THE CASE 1111   + + - -
!
      INTERACT=0
      IF(JA.EQ.JB) THEN
        IF(JJA.NE.JJB) THEN
          CALL RECOONESCALAR(NS,JA,JA,JA,JA,0,IAT)
          IF(IAT.EQ.0)RETURN
        END IF
!GG INT
        INTERACT = 1
        RETURN
!
!        CALL PERKO2(JA,JA,JA,JA,1)
!        QM1=HALF
!        QM2=-HALF
!        CALL WJ1(IK1,BK1,ID1,BD1,0,QM1,QM2,WJ)
!        IF(DABS(WJ).GT.EPS) THEN
!           RECOUPL=ONE/DSQRT(DBLE(IK1(6)+1))
!           COEFF=WJ*RECOUPL*DSQRT(DBLE(ID1(3)+1))
!           COEFF=-COEFF
!        END IF
      END IF
      RETURN
      END SUBROUTINE ONESCALAR1INT
