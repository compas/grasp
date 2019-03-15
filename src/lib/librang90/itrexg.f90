!*******************************************************************
!                                                                  *
      INTEGER FUNCTION ITREXG(I1,I2,I3,I4,K)
!
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: I1, I2, I3, I4
      INTEGER, INTENT(OUT) :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J
!-----------------------------------------------
      J=MAX0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG=J
      RETURN
      END FUNCTION ITREXG
