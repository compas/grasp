      MODULE hblock_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:29:39  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NBLOCK  !! NCFBLK as an integer needs to be renamed
!     REAL(DOUBLE) :: PNCFBLK
!     REAL(DOUBLE) :: PNEVBLK, PNCMAXBLK
      INTEGER, DIMENSION(:), pointer :: ncfblk  !! this is a problem
      INTEGER, DIMENSION(:), pointer :: nevblk, ncmaxblk
      END MODULE hblock_C
