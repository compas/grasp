      MODULE wave_C
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, DIMENSION(NNNW) :: MF
      REAL(DOUBLE), DIMENSION(NNNW) :: PZ
      INTEGER, DIMENSION(NNNW) :: MFFF
      REAL(DOUBLE), DIMENSION(NNNW) :: PZFF
      INTEGER, DIMENSION(NNNW) :: MFII
      REAL(DOUBLE), DIMENSION(NNNW) :: PZII
      REAL(DOUBLE), DIMENSION(:,:), pointer :: PF,QF
      REAL(DOUBLE), DIMENSION(:,:), pointer :: pfff,qfff,pfii,qfii
      END MODULE wave_C
