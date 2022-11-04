      MODULE natorbnew_I   
      INTERFACE
!...Translated by Gediminas Gaigalas 11/18/19
      SUBROUTINE natorbnew(NAME,DENS1VEC,DINT1VEC)
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW, NNNP
      USE prnt_C,          ONLY : NVEC
      REAL(DOUBLE), DIMENSION(NVEC,NNNP), INTENT(OUT) :: DENS1VEC
      CHARACTER*24, INTENT(IN) :: NAME
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NNNP), INTENT(IN) :: DINT1VEC
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
