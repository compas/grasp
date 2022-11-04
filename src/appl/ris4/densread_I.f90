      MODULE densread_I
      INTERFACE
!...Translated by Gediminas Gaigalas 11/18/19
      SUBROUTINE densread (DINT1,DINT2,DINT3,                          &
                           DINT4,DINT5,DINT6,                          &
                           DINT7)
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW, NNNP
      USE prnt_C,           ONLY: NVEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      END SUBROUTINE
      END INTERFACE
      END MODULE
