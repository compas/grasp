!*******************************************************************
!                                                                  *
      INTEGER FUNCTION JTHN(K,N,I)
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, N, K
!-----------------------------------------------
      IF(N == 1) THEN
        JTHN=MOD(K,I)
      ELSEIF(N == 2) THEN
        JTHN=MOD(K/I,I)
      ELSEIF(N == 3) THEN
        JTHN=MOD(K/(I*I),I)
      ELSEIF(N == 4) THEN
        JTHN=MOD(K/(I*I*I),I)
      ELSEIF(N == 5) THEN
        JTHN=MOD(K/(I*I*I*I),I)
      ELSE
        WRITE(6,'(A)') ' ERROR IN JTHN '
        STOP
      ENDIF
      RETURN
      END FUNCTION JTHN
