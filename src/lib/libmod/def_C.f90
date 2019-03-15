!
!***********************************************************************
!                                                                      *
      MODULE def_C
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNP
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE) :: TENMAX, EXPMAX, EXPMIN, PRECIS
      REAL(DOUBLE) :: AUCM, AUEV, CCMS, FASI, FBSI
      REAL(DOUBLE) :: FMTOAU, AUMAMU, B1
      INTEGER :: IONCTY, NELEC
      REAL(DOUBLE) :: EMN, Z
      INTEGER :: IONCTYFF, NELECFF
      REAL(DOUBLE) :: EMNFF, ZFF
      INTEGER :: IONCTYII, NELECII
      REAL(DOUBLE) :: EMNII, ZII
      INTEGER :: NELECR
      REAL(DOUBLE) :: C
      REAL(DOUBLE) :: EMPAM, RBCM
      INTEGER :: NSCF, NSIC, NSOLV
      REAL(DOUBLE) :: ACCY
      REAL(DOUBLE), DIMENSION(:), pointer :: wt, weight
      INTEGER :: NCMIN, NCMAX
      INTEGER, DIMENSION(:), pointer :: iccmin
      REAL(DOUBLE) :: CVAC, PI
      REAL(DOUBLE), DIMENSION(NNNP) :: DP, DQ
      END MODULE def_C
