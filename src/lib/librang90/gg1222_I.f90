      MODULE GG1222_I
      INTERFACE
!
      SUBROUTINE GG1222(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,K1, &
                        QM1,QM2,QM3,QM4,WW)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK1, IK2, ID1, ID2
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK1, BK2, BD1, BD2
      REAL(DOUBLE), INTENT(OUT)              :: WW
      END SUBROUTINE
      END INTERFACE
      END MODULE
