      MODULE scf_C 
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNW
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE), DIMENSION(NNNW) :: UCF 
      INTEGER, DIMENSION(NNNW) :: METHOD 
      REAL(DOUBLE), DIMENSION(NNNW) :: SCNSTY 
      REAL(DOUBLE) :: EPSMIN, EPSMAX, EMIN, EMAX, ZINF 
      INTEGER :: NDCOF, NXCOF, NYCOF, NDDIM, NXDIM, NYDIM 
      REAL(DOUBLE), DIMENSION(:), pointer :: da, xa, ya 
      INTEGER, DIMENSION(:), pointer :: nda, nxa, nya
      END MODULE scf_C 
