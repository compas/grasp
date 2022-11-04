      MODULE printaLS_I
      INTERFACE
      SUBROUTINE PRINTALS (INUM_II,INUM_FF,ASFA,ASFB,I,J,OMEGA,FACTOR)
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER  :: I, J, INUM_II, INUM_FF
      REAL(DOUBLE), INTENT(IN) :: ASFA, ASFB
      REAL(DOUBLE), INTENT(IN) :: OMEGA, FACTOR
      END SUBROUTINE
      END INTERFACE
      END MODULE
