!*******************************************************************
!                                                                  *
      SUBROUTINE EL4INT(JJA,JJB,JA,JB,JC,JD,ICOLBREI,INTERACT)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!     SUBROUTINE CALLED: EL41                                      *
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
      USE el41INT_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
      INTEGER, INTENT(OUT) :: INTERACT
!-----------------------------------------------
      INTERACT=0
      IF(NPEEL.LE.2)RETURN
      IF(JA.EQ.JB) THEN
        CALL EL41INT(JJA,JJB,JC,JD,JA,1,JA,JB,JC,JD,ICOLBREI,INTERACT)
      ELSE IF(JC.EQ.JD) THEN
        CALL EL41INT(JJA,JJB,JA,JB,JC,2,JA,JB,JC,JD,ICOLBREI,INTERACT)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL4INT ')
      END SUBROUTINE EL4INT
