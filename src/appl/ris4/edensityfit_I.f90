      MODULE edensityfit_I
      INTERFACE
!...Translated by Gediminas Gaigalas 11/18/19
      SUBROUTINE edensityfit (XVEC,YVEC,Z,PAR,NRNUC,F,RHO,RES)
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP, NNN1
      REAL(DOUBLE), DIMENSION(NNN1) :: XVEC
      REAL(DOUBLE), DIMENSION(NNNP) :: YVEC
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      REAL(DOUBLE), DIMENSION(5) :: P, F
      REAL(DOUBLE), DIMENSION(2) :: PAR
      REAL(DOUBLE) :: Z, DRMS, RHO, RES
      INTEGER      :: NRNUC, NPARFIT
      END SUBROUTINE
      END INTERFACE
      END MODULE
