!
!***********************************************************************
!                                                                      *
      MODULE orb_C
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
      USE parameter_def,   ONLY:  NNNW
!...Created by Pacific-Sierra Research 77to90  4.3E  08:57:22  12/25/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17

      CHARACTER(LEN=2), DIMENSION(NNNW) :: NH
      CHARACTER(LEN=2), DIMENSION(NNNW) :: NHR
      REAL(DOUBLE), DIMENSION(NNNW) :: E, GAMA, PED
      INTEGER :: NCF, NW, NCFR, NWR
      INTEGER(BYTE), DIMENSION(:,:), pointer :: IQA
      INTEGER, DIMENSION(:,:), pointer :: IQAR
      INTEGER, DIMENSION(NNNW) :: NP, NAK, NPR, NAKR
      INTEGER, DIMENSION(NNNW) :: NKL, NKJ
      END MODULE orb_C
