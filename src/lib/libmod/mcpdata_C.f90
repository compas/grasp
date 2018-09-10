      MODULE mcpdata_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NCOEFF, NINTG
      INTEGER, DIMENSION(:), pointer :: jann, jbnn, intgrl, intptr
      REAL(DOUBLE), DIMENSION(:), pointer :: cnn
      END MODULE mcpdata_C 
