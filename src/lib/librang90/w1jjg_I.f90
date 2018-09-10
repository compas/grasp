      MODULE W1JJG_I
      INTERFACE
!
      SUBROUTINE W1JJG(K1,QM1,QM2,IK,BK,ID,BD,WW)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WW
      END SUBROUTINE
      END INTERFACE
      END MODULE
