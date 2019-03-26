      MODULE vpint_I
      INTERFACE
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      SUBROUTINE VPINT (IA,IB,TEGRAL)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER , INTENT(INOUT) :: IA, IB
      REAl(DOUBLE), INTENT(OUT) :: TEGRAL
      END SUBROUTINE vpint
      END INTERFACE
      END MODULE
