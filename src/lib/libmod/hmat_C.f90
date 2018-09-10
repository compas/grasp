      MODULE hmat_C 
      USE vast_kind_param, ONLY:  DOUBLE, LONG
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER(LONG) :: NELMNT 
      REAL(DOUBLE), DIMENSION(:), pointer :: emt
      INTEGER, DIMENSION(:), pointer :: iendc   
      INTEGER, DIMENSION(:), pointer :: irow
      INTEGER, DIMENSION(6) :: NTPITMP 
      INTEGER(LONG) :: NELMNTTMP
      INTEGER :: NCOEITMP, NCOECTMP, NCTEITMP, NCTECTMP, NMCBPTMP,  &
                 NCORETMP, NVPITMP, NKEITMP, NVINTITMP, NCFTMP 
      REAL(DOUBLE) :: CUTOFFTMP 
      END MODULE hmat_C 
