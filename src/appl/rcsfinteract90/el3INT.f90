!*******************************************************************
!                                                                  *
      SUBROUTINE EL3INT(JJA,JJB,JA,JB,JC,JD,ICOLBREI,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     SUBROUTINE CALLED: EL31,EL32,EL33                            *
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
      USE el31INT_I
      USE el32INT_I
      USE el33INT_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
      IF(NPEEL.LE.1)RETURN
      IF(JB.EQ.JD) THEN
        IF(JA.EQ.JB.OR.JC.EQ.JB) THEN
          IF(JA.EQ.JC)GO TO 10
          IF(JC.NE.JB) THEN
            CALL EL32INT(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI,INTERACT)
          ELSE
            CALL EL31INT(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI,INTERACT)
          END IF
        ELSE
          CALL EL33INT(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI,INTERACT)
        END IF
        RETURN
      ELSE IF(JA.EQ.JC) THEN
        IF(JB.EQ.JA.OR.JD.EQ.JA) THEN
          IF(JB.EQ.JD)GO TO 10
          IF(JD.NE.JA) THEN
            CALL EL32INT(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI,INTERACT)
          ELSE
            CALL EL31INT(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI,INTERACT)
          END IF
        ELSE
          CALL EL33INT(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI,INTERACT)
        END IF
        RETURN
      ELSE IF(JA.EQ.JD) THEN
        IF(JB.EQ.JA.OR.JC.EQ.JA) THEN
          IF(JB.EQ.JC)GO TO 10
          IF(JC.NE.JD) THEN
             CALL EL32INT(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI,INTERACT)
          ELSE
             CALL EL31INT(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI,INTERACT)
          END IF
        ELSE
          CALL EL33INT(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI,INTERACT)
        END IF
        RETURN
      ELSE IF(JB.EQ.JC) THEN
        IF(JA.EQ.JB.OR.JD.EQ.JB) THEN
          IF(JA.EQ.JD)GO TO 10
          IF(JD.NE.JB) THEN
            CALL EL32INT(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI,INTERACT)
          ELSE
            CALL EL31INT(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI,INTERACT)
          END IF
        ELSE
          CALL EL33INT(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI,INTERACT)
        END IF
        RETURN
      END IF
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL3INT  PMGG RAGG')
      STOP
      END SUBROUTINE EL3INT
