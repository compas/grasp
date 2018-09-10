!*******************************************************************
!                                                                  *
      SUBROUTINE EL5INT(JJA,JJB,JA,JB,JC,JD,ICOLBREI,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!      SUBROUTINE CALLED: EL51INT,EL52INT,EL53INT                  *
!                                                                  *
!   Written by  G. Gaigalas                   NIST, December 2015  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE m_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE el51INT_I
      USE el52INT_I
      USE el53INT_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
      IF(NPEEL.LE.3)RETURN
      IF(JB.LT.JC) THEN
        CALL EL51INT(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI,INTERACT)
      ELSE IF(JA.GT.JD.AND.JB.GT.JD) THEN
        CALL EL51INT(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI,INTERACT)
      ELSE IF(JB.GT.JC.AND.JB.LT.JD.AND.JA.LT.JC) THEN
        CALL EL52INT(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI,INTERACT)
      ELSE IF(JB.GT.JC.AND.JB.GT.JD.AND.JA.GT.JC) THEN
        CALL EL52INT(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI,INTERACT)
      ELSE IF(JB.GT.JC.AND.JB.GT.JD.AND.JA.LT.JC) THEN
        CALL EL53INT(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI,INTERACT)
      ELSE IF(JB.GT.JC.AND.JB.LT.JD.AND.JA.GT.JC) THEN
        CALL EL53INT(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI,INTERACT)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL5INT ')
      END SUBROUTINE EL5INT
