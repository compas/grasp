!*******************************************************************
!                                                                  *
      SUBROUTINE RMEW7JJ(J1,J2,K1,K2,COEF)
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rmew7bjj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(DOUBLE), INTENT(OUT) :: COEF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(14) :: I00, I10, I01
!-----------------------------------------------
      DATA I00/32,48,128,80,96,128,20,60,20,108,36,44,156,68/
      DATA I10/6,9,120,15,18,24,2*30,0,54,2*0,78,0/
      DATA I01/10,35,168,165,286,680,0,30,10,180,60,110,546,408/
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 12 .OR. J1 > 25) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I00(J1-11)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I10(J1-11)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I01(J1-11))/DBLE(7))
      ELSE
        CALL RMEW7BJJ(J1,J2,K1,K2,COEF)
      ENDIF
      RETURN
      END SUBROUTINE RMEW7JJ
