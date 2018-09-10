!*******************************************************************
!                                                                  *
      SUBROUTINE A1JJ(IK,BK,ID,BD,QM1,A)
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
      USE CONS_C,          ONLY: ZERO, EPS
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(IN)               :: QM1
      REAL(DOUBLE), INTENT(OUT)              :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IFAZ, ISUMA
      REAL(DOUBLE) :: AB
!-----------------------------------------------
      A=ZERO
      IF(QM1 < EPS) THEN
        ISUMA=(ID(6)+1)*ID(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IFAZ=ID(6)+ID(3)-IK(6)+ID(4)*2
        IF((IFAZ/4)*4 /= IFAZ)A=-A
      ELSE
        ISUMA=(IK(6)+1)*IK(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IF((IK(4)/2)*2 /= IK(4))A=-A
      ENDIF
      RETURN
      END SUBROUTINE A1JJ
