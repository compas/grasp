      MODULE iabint_I
      INTERFACE
      SUBROUTINE IABINT (IA,IB,TEGRAL)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(INOUT) :: ia, ib
      REAL(DOUBLE), INTENT(out) :: tegral
      END SUBROUTINE
      END INTERFACE
      END MODULE
