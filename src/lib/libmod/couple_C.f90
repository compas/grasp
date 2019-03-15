      MODULE couple_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, PARAMETER :: MANGM = 60
      INTEGER, PARAMETER :: MTRIAD = 12
      INTEGER, DIMENSION(MANGM) :: J1
      INTEGER, DIMENSION(MTRIAD,3) :: J2, J3
      INTEGER :: MJA, NJA
      LOGICAL, DIMENSION(MANGM) :: FREE
      END MODULE couple_C
