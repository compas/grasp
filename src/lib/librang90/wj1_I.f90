      MODULE WJ1_I
      INTERFACE
!
      SUBROUTINE WJ1(IK,BK,ID,BD,K2,QM1,QM2,WJ)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K2
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WJ
      END SUBROUTINE
      END INTERFACE
      END MODULE
