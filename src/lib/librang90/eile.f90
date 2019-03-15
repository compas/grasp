!*******************************************************************
!                                                                  *
      SUBROUTINE EILE(JA,JB,JC,JAA,JBB,JCC)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 02  ------------   *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA, JB, JC
      INTEGER, INTENT(OUT) :: JAA, JBB, JCC
!-----------------------------------------------
      JAA=JA
      JCC=JA
      IF(JAA > JB)JAA=JB
      IF(JCC < JB)JCC=JB
      IF(JAA > JC)JAA=JC
      IF(JCC < JC)JCC=JC
      IF((JA > JAA).AND.(JA < JCC))JBB=JA
      IF((JB > JAA).AND.(JB < JCC))JBB=JB
      IF((JC > JAA).AND.(JC < JCC))JBB=JC
      RETURN
      END SUBROUTINE EILE
