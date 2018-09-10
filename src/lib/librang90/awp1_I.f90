      MODULE AWP1_I
      INTERFACE
!
      SUBROUTINE AWP1( IK, BK, ID, BD, K1, BK2, QM1, QM2, QM3, AW)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: BK2, QM1, QM2, QM3
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: AW
      END SUBROUTINE
      END INTERFACE
      END MODULE
