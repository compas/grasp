      MODULE sbdat1_C
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, PARAMETER :: NLMAX = 20
      INTEGER, DIMENSION(NLMAX,NLMAX) :: NSHLP
      INTEGER, DIMENSION(NLMAX,NNNW) :: NSHLPP
      END MODULE sbdat1_C
