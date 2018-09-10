      MODULE A1JJ_I
      INTERFACE
!
      SUBROUTINE A1JJ (IK, BK, ID, BD, QM1, A)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(IN)               :: QM1
      REAL(DOUBLE), INTENT(OUT)              :: A
      END SUBROUTINE
      END INTERFACE
      END MODULE
