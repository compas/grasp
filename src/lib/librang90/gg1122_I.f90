      MODULE GG1122_I
      INTERFACE
!
      SUBROUTINE GG1122(K1,K2,QM1,QM2,QM3,QM4,AA)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)               :: K1, K2
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(OUT)              :: AA
      END SUBROUTINE
      END INTERFACE
      END MODULE
