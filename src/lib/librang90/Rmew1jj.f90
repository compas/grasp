!*******************************************************************
!                                                                  *
      SUBROUTINE RMEW1JJ(J1,J2,K1,K2,COEF)
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
      USE CONS_C,          ONLY: ZERO
      USE ribojj_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(DOUBLE), INTENT(OUT) :: COEF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(2) :: I10, I01
!-----------------------------------------------
      DATA I01/6,0/
      DATA I10/0,6/
!
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 > 2) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
        COEF=-DSQRT(DBLE(2))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
        COEF=-DSQRT(DBLE(I10(J1)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
        COEF=-DSQRT(DBLE(I01(J1)))
      ELSE
        WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
        WRITE(0,'(A)') ' ERROR IN SUB. RMEW1JJ '
        STOP
      ENDIF
      RETURN
      END SUBROUTINE RMEW1JJ
