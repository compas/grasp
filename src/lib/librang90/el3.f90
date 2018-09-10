!*******************************************************************
!                                                                  *
      SUBROUTINE EL3(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     SUBROUTINE CALLED: EL31,EL32,EL33                            *
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
      USE m_C,         ONLY:  NPEEL
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE el31_I
      USE el32_I
      USE el33_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 1)RETURN
      IF(JB == JD) THEN
        IF(JA == JB.OR.JC == JB) THEN
          IF(JA == JC)GO TO 10
          IF(JC /= JB) THEN
            CALL EL32(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JC) THEN
        IF(JB == JA.OR.JD == JA) THEN
          IF(JB == JD)GO TO 10
          IF(JD /= JA) THEN
            CALL EL32(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JA == JD) THEN
        IF(JB == JA.OR.JC == JA) THEN
          IF(JB == JC)GO TO 10
          IF(JC /= JD) THEN
             CALL EL32(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          ELSE
             CALL EL31(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      ELSE IF(JB == JC) THEN
        IF(JA == JB.OR.JD == JB) THEN
          IF(JA == JD)GO TO 10
          IF(JD /= JB) THEN
            CALL EL32(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          ELSE
            CALL EL31(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
          END IF
        ELSE
          CALL EL33(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI)
        END IF
        RETURN
      END IF
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL3  PMGG RAGG')
      STOP
      END SUBROUTINE EL3
