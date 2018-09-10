      MODULE vpilst_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:35:13  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NDVPA, NVPI 
      INTEGER, DIMENSION(:), pointer :: indvpi
      REAL(DOUBLE), DIMENSION(:), pointer :: valvpi
      LOGICAL :: FRSTVP 
      END MODULE vpilst_C 
