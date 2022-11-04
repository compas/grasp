      MODULE sigma_1_I
      INTERFACE
      SUBROUTINE sigma_1 (IPAR,I1,I2,APART)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN) :: IPAR, I1, I2
      REAL(DOUBLE), INTENT(OUT) :: APART
      END SUBROUTINE sigma_1
      END INTERFACE
      END MODULE
