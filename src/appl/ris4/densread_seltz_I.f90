      MODULE densread_seltz_I
      INTERFACE
!...Translated by Gediminas Gaigalas 11/18/19
      SUBROUTINE densread_seltz (DINT1,DINT2,DINT3,                    &
                           DINT4,DINT5,DINT6,                          &
                           DINT7,DINT1VEC,DENS1VEC,NRNUC)
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW, NNNP
      USE prnt_C,           ONLY: NVEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      INTEGER, INTENT(IN) :: NRNUC
      REAL(DOUBLE), DIMENSION(NVEC,NRNUC), INTENT(OUT)     :: DENS1VEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NRNUC), INTENT(IN) :: DINT1VEC
      END SUBROUTINE
      END INTERFACE
      END MODULE
