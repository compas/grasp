!***********************************************************************
!                                                                      *
      SUBROUTINE MOHR(N, KAPPA, Z, FZALFA)
!                                                                      *
!   The  function  F (Z*alpha)  for the  1s  2s  2p-  2p  symmetries   *
!   is computed here.    A value is obtained by interpolating in, or   *
!   extrapolating from, the table due to  P J Mohr.   See  P J Mohr,   *
!   At Data Nucl Data Tables 29 (1983) 453.                            *
!                                                                      *
!   Call(s) to: [LIB92]: INTERP.                                       *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE interp_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: KAPPA
      REAL(DOUBLE)  :: Z
      REAL(DOUBLE) , INTENT(OUT) :: FZALFA
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!   Number of data points
      INTEGER, PARAMETER :: NUMVAL = 12
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NUMVAL) :: VAL1S, VAL2S, VAL2P1, VAL2P3, ARG
      REAL(DOUBLE) :: ACCY, VALUE
!-----------------------------------------------
!
!
!   1s data :
!
      DATA VAL1S/ 10.3168D00, 4.6540D00, 3.2460D00, 2.5519D00, 2.1351D00, &
         1.8644D00, 1.6838D00, 1.5675D00, 1.5032D00, 1.4880D00, 1.5317D00, &
         1.6614D00/
!
!   2s data:
!
      DATA VAL2S/ 10.5468D00, 4.8930D00, 3.5063D00, 2.8391D00, 2.4550D00, &
         2.2244D00, 2.0948D00, 2.0435D00, 2.0650D00, 2.1690D00, 2.3870D00, &
         2.7980D00/
!
!   2p- data:
!
      DATA VAL2P1/ -0.1264D00, -0.1145D00, -0.0922D00, -0.0641D00, -0.0308D00, &
         0.0082D00, 0.0549D00, 0.1129D00, 0.1884D00, 0.2934D00, 0.4530D00, &
         0.7250D00/
!
!   2p data:
!
      DATA VAL2P3/ 0.1235D00, 0.1303D00, 0.1436D00, 0.1604D00, 0.1794D00, &
         0.1999D00, 0.2215D00, 0.2440D00, 0.2671D00, 0.2906D00, 0.3141D00, &
         0.3367D00/
!
!   Z data:
!
      DATA ARG/ 1.0D00, 10.0D00, 20.0D00, 30.0D00, 40.0D00, 50.0D00, 60.0D00, &
         70.0D00, 80.0D00, 90.0D00, 100.0D00, 110.0D00/
!
!----------------------------------------------------------------------*
!
!   Convergence criterion for interpolation
!
      DATA ACCY/ 1.0D-03/
!
!   Interpolate or issue error message as appropriate
!
      IF (N == 1) THEN
         IF (KAPPA == (-1)) THEN
            CALL INTERP (ARG, VAL1S, NUMVAL, Z, VALUE, ACCY)
         ELSE
            WRITE (*, 300)
            WRITE (*, 301) N, KAPPA
            STOP
         ENDIF
      ELSE IF (N == 2) THEN
         SELECT CASE (KAPPA)
         CASE (-1)
            CALL INTERP (ARG, VAL2S, NUMVAL, Z, VALUE, ACCY)
         CASE (1)
            CALL INTERP (ARG, VAL2P1, NUMVAL, Z, VALUE, ACCY)
         CASE (-2)
            CALL INTERP (ARG, VAL2P3, NUMVAL, Z, VALUE, ACCY)
         CASE DEFAULT
            WRITE (*, 300)
            WRITE (*, 301) N, KAPPA
            STOP
         END SELECT
      ELSE
         WRITE (*, 300)
         WRITE (*, 302) N
         STOP
      ENDIF
!
      FZALFA = VALUE
!
      RETURN
!
  300 FORMAT('MOHR:')
  301 FORMAT(' Principal quantum number, ',I12,', kappa, ',1I3,'.')
  302 FORMAT(' Principal quantum number, ',1I2,', Should be either 1 or 2.')
      RETURN
!
      END SUBROUTINE MOHR
