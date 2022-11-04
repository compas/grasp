      MODULE smsnew_I
      INTERFACE
!...Translated by Gediminas Gaigalas 11/18/19
      SUBROUTINE smsnew (DOIT,VINT,VINT2)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      INTEGER, INTENT(IN) :: DOIT
      END SUBROUTINE
      END INTERFACE
      END MODULE
