      MODULE int_C 
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNP
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: MTP0 
      REAL(DOUBLE), DIMENSION(NNNP) :: P, Q 
      REAL(DOUBLE) :: P0, Q0 
      REAL(DOUBLE), DIMENSION(NNNP) :: TF, TG, XU, XV 
      END MODULE int_C 
