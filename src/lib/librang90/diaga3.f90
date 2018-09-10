!*******************************************************************
!                                                                  *
      SUBROUTINE DIAGA3(JA1,JA2,KA,IRE,IAT,REC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 03  --------------  *
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
      USE CONS_C,          ONLY: ZERO, ONE
      USE m_C,             ONLY: JLIST, JJQ1, JJC1, JJC2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA1, JA2, KA, IRE
      INTEGER, INTENT(OUT) :: IAT
      REAL(DOUBLE), INTENT(OUT) :: REC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,LL1,JI,KK1,KK2,ITI,ITI1,ITIS,ITI1S,IFAZ
      REAL(DOUBLE) :: AA, A3
!-----------------------------------------------
      REC = ZERO
      AA=ONE
      I=JA1+1
      IF(JA1 == 1)I=I+1
      IF(I >= JA2) THEN
        REC=AA
        IAT=1
      ELSE
    1   LL1=JLIST(I)
        JI=JJQ1(3,LL1)-1
        KK2=I-2
        ITI=JJC1(KK2)-1
        ITIS=JJC2(KK2)-1
        KK1=I-1
        ITI1=JJC1(KK1)-1
        ITI1S=JJC2(KK1)-1
        IF(IRE /= 0) THEN
          CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
          A3=A3*SQRT(DBLE((ITI+1)*(ITI1S+1)))
          IFAZ=KA+JI+ITI+ITI1S
          IF((IFAZ/4)*4 /= IFAZ)A3=-A3
          AA=AA*A3
        ELSE
          IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) == 0)RETURN
        END IF
        I=I+1
        IF(I == JA2) THEN
          REC=AA
          IAT=1
        ELSE
          GO TO 1
        END IF
      END IF
      RETURN
      END SUBROUTINE DIAGA3
