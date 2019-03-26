!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION SKFUN (K, X)
!                                                                      *
!   Computes the function                                              *
!                                            n  nx                     *
!                             infinity   (-1)  e                       *
!                     S (x) =   Sum      -------                       *
!                               n=1          k                         *
!                                           n                          *
!                                                                      *
!   See, for instance, F. A. Parpia and A. K. Mohanty ``Relativistic   *
!   basis set calculations for atoms with Fermi nuclei'', Phys Rev A   *
!   (1992) in press.                                                   *
!                                                                      *
!   Call(s) to: SKFUN.                                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:31:28   1/ 3/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K
      REAL(DOUBLE) , INTENT(IN) :: X
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: QUASIZERO = 1.D-15
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: DNUMER, EN, BASE, DELTA
!-----------------------------------------------
!xhe

      DNUMER = 1.0D00
      EN = 0.0D00
      BASE = -EXP(X)
      SKFUN = 0.0D00
      DNUMER = DNUMER*BASE
      EN = EN + 1.0D00
      DELTA = DNUMER/EN**K
      SKFUN = SKFUN + DELTA
      DO WHILE(ABS(DELTA/SKFUN) > QUASIZERO)
         DNUMER = DNUMER*BASE
         EN = EN + 1.0D00
         DELTA = DNUMER/EN**K
         SKFUN = SKFUN + DELTA
!xhb
!xh      IF (ABS (DELTA/SKFUN) .GT. 1.0D-15 ) GOTO 1
      END DO
!xhe
!
      RETURN
      END FUNCTION SKFUN
