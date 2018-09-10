      MODULE brint1_I   
      INTERFACE
!...Translated by Charlotte Froese Fischer 
!                       Gediminas Gaigalas  10/05/17
      SUBROUTINE BRINT1 (IA,IB,IC,ID,K,TEGRAL)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: ia, ib, ic, id, k
      REAL(DOUBLE), INTENT(out) :: tegral
      END SUBROUTINE
      END INTERFACE
      END MODULE 
