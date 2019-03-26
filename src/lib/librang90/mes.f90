!*******************************************************************
!                                                                  *
      SUBROUTINE MES(I)
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J
      CHARACTER(LEN=10), DIMENSION(6) :: STRING5
!
      DATA STRING5/'   I T L S',' I T L S 2',' I T L S 3', &
                   'AW P 1 L S','WA P 1 L S','       W 1'/
!-----------------------------------------------
      IF(I > 50) THEN
        J=I-50
        WRITE(6,'(A)') ' error in func./sub. '
        WRITE(6,'(20X,A10)') STRING5(J)
        WRITE(6,'(A)') ' susimaise f sluoksnio termu kodavimas  '
      ELSE
        WRITE(6,'(A)') ' yra daugiau nei 2 ele. sluoks. f,g,h,i,k,l,m'
        WRITE(6,'(3X,I5)') I
        IF(I == 1) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 G '
        ELSEIF(I == 2) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 '
        ELSEIF(I == 11) THEN
          WRITE(6,'(A)') ' error in Function     I T L S  '
        ELSEIF(I == 12) THEN
          WRITE(6,'(A)') ' error in Function     I T L S 2  '
        ELSEIF(I == 13) THEN
          WRITE(6,'(A)') ' error in Function     I T L S 3  '
        ELSEIF(I == 30) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 A 3 A 4 L S '
        ELSEIF(I == 31) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 L S '
        ELSEIF(I == 32) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
        ELSEIF(I == 33) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
        ELSEIF(I == 34) THEN
          WRITE(6,'(A)') ' error in Subroutine   W A 1 A 2 L S '
        ELSEIF(I == 35) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 W 2 L S '
        ELSE
          WRITE(6,'(A)') ' error in unknown Subroutine  '
        ENDIF
      ENDIF
      WRITE(6,'(A)') ' Contact to   G. Gediminas please ! ! ! '
      STOP
      END SUBROUTINE MES
