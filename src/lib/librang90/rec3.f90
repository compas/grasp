!*******************************************************************
!                                                                  *
      SUBROUTINE REC3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 07  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  RECO3                                    *
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reco3_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: JA1,JA2,JA3,K1,K2,K3,IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IFAZ
!-----------------------------------------------
      IF((JA3 > JA1).AND.(JA3 > JA2)) THEN
        IF(JA1-JA2 < 0) THEN
          CALL RECO3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,RECC)
        ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA2,JA1,JA3,K2,K1,K3,IRE,IAT,RECC)
          IFAZ=K1+K2-K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
        ELSE
          GO TO 10
        END IF
      ELSE IF((JA3 < JA1).AND.(JA3 < JA2)) THEN
        IF(JA1-JA2 < 0) THEN
          CALL RECO3(JA3,JA1,JA2,K3,K1,K2,IRE,IAT,RECC)
          IF((K3/2)*2 /= K3)RECC=-RECC
        ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA3,JA2,JA1,K3,K2,K1,IRE,IAT,RECC)
          IFAZ=K1+K2+K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
        ELSE
          GO TO 10
        END IF
      ELSE
        IF(JA1-JA2 < 0)THEN
          CALL RECO3(JA1,JA3,JA2,K1,K3,K2,IRE,IAT,RECC)
          IFAZ=K1-K2-K3
          IF((IFAZ/4)*4 /= IFAZ)RECC=-RECC
        ELSE IF(JA1-JA2 > 0) THEN
          CALL RECO3(JA2,JA3,JA1,K2,K3,K1,IRE,IAT,RECC)
          IF((K1/2)*2 /= K1)RECC=-RECC
        ELSE
          GO TO 10
        END IF
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN REC')
      STOP
      END SUBROUTINE REC3
