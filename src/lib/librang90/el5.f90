!*******************************************************************
!                                                                  *
      SUBROUTINE EL5(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!      SUBROUTINE CALLED: EL51,EL52,EL53                           *
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
      USE m_C,          ONLY: NPEEL
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE el51_I
      USE el52_I
      USE el53_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JJA,JJB,JA,JB,JC,JD,ICOLBREI
!-----------------------------------------------
      IF(NPEEL <= 3)RETURN
      IF(JB < JC) THEN
        CALL EL51(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI)
      ELSE IF(JA > JD.AND.JB > JD) THEN
        CALL EL51(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA < JC) THEN
        CALL EL52(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA > JC) THEN
        CALL EL52(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI)
      ELSE IF(JB > JC.AND.JB > JD.AND.JA < JC) THEN
        CALL EL53(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI)
      ELSE IF(JB > JC.AND.JB < JD.AND.JA > JC) THEN
        CALL EL53(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI)
      ELSE
        WRITE(99,100)
        STOP
      END IF
      RETURN
  100 FORMAT(5X,'ERRO IN EL5 ')
      END SUBROUTINE EL5
