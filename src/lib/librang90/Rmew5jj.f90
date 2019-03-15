!*******************************************************************
!                                                                  *
      SUBROUTINE RMEW5JJ(J1,J2,K1,K2,COEF)
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
      INTEGER                 :: JI1, JI2
      INTEGER, DIMENSION(6)   :: I00, I01
      INTEGER, DIMENSION(3,3) :: I12P,I12A,I03P, &
                       I03A,I14P,I14A,I05P,I05A
!-----------------------------------------------
      DATA I00/54,12,30,12,30,54/
      DATA I01/126,12,198,0,48,288/
      DATA I12P/420,360,270,360,2*0,-270,2*0/
      DATA I12A/0,1960,0,-1960,-1000,2430,0,2430,1980/
      DATA I03P/-882,3*0,384,-400,0,400,286/
      DATA I03A/4*0,162,-300,0,-300,22/
      DATA I14P/3780,-720,-4950,-720,2*0,4950,2*0/
      DATA I14A/2*0,3528,0,2430,1980,-3528,1980,-7722/
      DATA I05P/-1386,4*0,-440,0,440,-1430/
      DATA I05A/5*0,330,0,330,572/
!
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 6 .OR. J1 > 11) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I00(J1-5)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
        IF(J1 /= J2) RETURN
        IF(J1 == 6) THEN
           COEF=-DSQRT(DBLE(48))
        ELSEIF(J1 == 9) THEN
           COEF=-DSQRT(DBLE(20))
        ELSEIF(J1 == 10) THEN
           COEF=-DSQRT(DBLE(10))
        ELSEIF(J1 == 11) THEN
           COEF=-DSQRT(DBLE(18))
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
        IF(J1 /= J2) RETURN
        COEF=-DSQRT(DBLE(I01(J1-5))/DBLE(7))
      ELSEIF(K1 == 1 .AND. K2 == 2) THEN
        IF(J1 < 9) THEN
           JI1=J1-5
           JI2=J2-5
           IF(I12P(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I12P(JI1,JI2))/DBLE(7))
           ELSE
              COEF=-DSQRT(-DBLE(I12P(JI1,JI2))/DBLE(7))
           ENDIF
        ELSE
           JI1=J1-8
           JI2=J2-8
           IF(I12A(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I12A(JI1,JI2))/DBLE(49))
           ELSE
              COEF=-DSQRT(-DBLE(I12A(JI1,JI2))/DBLE(49))
           ENDIF
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
        IF(J1 < 9) THEN
           JI1=J1-5
           JI2=J2-5
           IF(I03P(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I03P(JI1,JI2))/DBLE(21))
           ELSE
              COEF=-DSQRT(-DBLE(I03P(JI1,JI2))/DBLE(21))
           ENDIF
        ELSE
           JI1=J1-8
           JI2=J2-8
           IF(I03A(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I03A(JI1,JI2))/DBLE(7))
           ELSE
              COEF=-DSQRT(-DBLE(I03A(JI1,JI2))/DBLE(7))
           ENDIF
        ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 4) THEN
        IF(J1 < 9) THEN
           JI1=J1-5
           JI2=J2-5
           IF(I14P(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I14P(JI1,JI2))/DBLE(35))
           ELSE
              COEF=-DSQRT(-DBLE(I14P(JI1,JI2))/DBLE(35))
           ENDIF
        ELSE
           JI1=J1-8
           JI2=J2-8
           IF(I14A(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I14A(JI1,JI2))/DBLE(49))
           ELSE
              COEF=-DSQRT(-DBLE(I14A(JI1,JI2))/DBLE(49))
           ENDIF
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 5) THEN
        IF(J1 < 9) THEN
           JI1=J1-5
           JI2=J2-5
           IF(I05P(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I05P(JI1,JI2))/DBLE(21))
           ELSE
              COEF=-DSQRT(-DBLE(I05P(JI1,JI2))/DBLE(21))
           ENDIF
        ELSE
           JI1=J1-8
           JI2=J2-8
           IF(I05A(JI1,JI2) >= 0) THEN
              COEF=DSQRT(DBLE(I05A(JI1,JI2))/DBLE(7))
           ELSE
              COEF=-DSQRT(-DBLE(I05A(JI1,JI2))/DBLE(7))
           ENDIF
        ENDIF
      ELSE
         WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
         WRITE(0,'(A)') ' ERROR IN SUB. RMEW5JJ '
         STOP
      ENDIF
      RETURN
      END SUBROUTINE RMEW5JJ
