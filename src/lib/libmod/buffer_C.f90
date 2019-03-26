      MODULE buffer_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  13:03:28   1/25/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NBDIM, NVCOEF
      INTEGER, DIMENSION(:,:), pointer :: label
      REAL(DOUBLE), DIMENSION(:), pointer :: coeff
      END MODULE buffer_C
