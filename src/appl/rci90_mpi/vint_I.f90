      MODULE vint_I
      INTERFACE
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE VINT (IA,IB,TEGRAL)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER , INTENT(IN) :: IA, IB
      REAl(DOUBLE), INTENT(OUT) :: TEGRAL
      END SUBROUTINE vint
      END INTERFACE
      END MODULE
