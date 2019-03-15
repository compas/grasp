!*******************************************************************
!                                                                  *
      SUBROUTINE RMEW3JJ(J1,J2,K1,K2,COEF)
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
      INTEGER               :: JI1
      INTEGER, DIMENSION(3) :: I00, I10, I01
!-----------------------------------------------
      DATA I00/16,6,10/
      DATA I10/12,12,0/
      DATA I01/12,0,12/
!
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 2 .OR. J1 > 5) RETURN
      JI1=J1-2
      IF(K1 == 0 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I00(JI1)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I10(JI1)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I01(JI1)))
      ELSEIF(K1 == 1 .AND. K2 == 2) THEN
        IF(J1 == 3 .AND. J2 == 3) THEN
           COEF=DSQRT(DBLE(60))
        ELSEIF(J1 == 5 .AND. J2 == 4) THEN
           COEF=DSQRT(DBLE(30))
        ELSEIF(J1 == 4 .AND. J2 == 5) THEN
           COEF=-DSQRT(DBLE(30))
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
        IF(J1 == 3 .AND. J2 == 3) THEN
           COEF=-DSQRT(DBLE(28))
        ELSEIF(J1 == 5 .AND. J2 == 5) THEN
           COEF=DSQRT(DBLE(28))
        ENDIF
      ELSE
        WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
        WRITE(0,'(A)') ' ERROR IN SUB. RMEW3JJ '
        STOP
      ENDIF
      RETURN
      END SUBROUTINE RMEW3JJ
