      MODULE bilst_C
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, DIMENSION(6) :: NDTPA, NTPI
      INTEGER, DIMENSION(:), pointer :: indtp1, indtp2, indtp3, indtp4, &
                                        indtp5, indtp6
      REAL(DOUBLE), DIMENSION(:), pointer :: valtp1, valtp2, valtp3, valtp4, &
                                        valtp5, valtp6
      LOGICAL, DIMENSION(6) :: FIRST
      END MODULE bilst_C
