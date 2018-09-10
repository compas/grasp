      MODULE bess_C 
      USE vast_kind_param, ONLY:  DOUBLE 
      USE parameter_def,   ONLY:  NNNP
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE), DIMENSION(2) :: WIJ 
      REAL(DOUBLE), DIMENSION(2,2,NNNP) :: BESSJ, BESSN 
      REAL(DOUBLE), DIMENSION(NNNP,3) :: BJ 
      REAL(DOUBLE), DIMENSION(NNNP) :: TC, TD 
      END MODULE bess_C 
