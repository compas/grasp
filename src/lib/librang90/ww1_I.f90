      MODULE WW1_I
      INTERFACE
!
      SUBROUTINE WW1(IK,BK,ID,BD,K2,QM1,QM2,QM3,QM4,WW)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K2
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT)              :: WW
      END SUBROUTINE
      END INTERFACE
      END MODULE
