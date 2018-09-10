!*******************************************************************
!                                                                  *
      SUBROUTINE RWJJ(J,J1,J2,K1,K2,COEF)
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
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE ribojj_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rmew1jj_I
      USE rmew3jj_I
      USE rmew5jj_I
      USE rmew7jj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J, J1, J2, K1, K2
      REAL(DOUBLE), INTENT(OUT) :: COEF
!-----------------------------------------------
      IF(J == 1) THEN
         CALL RMEW1JJ(J1,J2,K1,K2,COEF)
      ELSEIF(J == 3) THEN
         CALL RMEW3JJ(J1,J2,K1,K2,COEF)
      ELSEIF(J == 5) THEN
         CALL RMEW5JJ(J1,J2,K1,K2,COEF)
      ELSEIF(J == 7) THEN
         CALL RMEW7JJ(J1,J2,K1,K2,COEF)
!GG        ELSEIF(J.EQ.9) THEN
!GG          CALL SUWJJ(K1,K2,J,J1,J2,COEF)
      ELSE
         WRITE(0,'(A,I5)') ' KLAIDA SUB. RWJJ J=',J
         STOP
      ENDIF
      RETURN
      END SUBROUTINE RWJJ
