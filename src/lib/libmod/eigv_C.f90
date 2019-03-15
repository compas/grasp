!
!***********************************************************************
!                                                                      *
      MODULE eigv_C
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE) :: EAV, EAVFF, EAVII
      REAL(DOUBLE), DIMENSION(:), pointer :: EVAL, EVALFF, EVALII
      REAL(DOUBLE), DIMENSION(:), pointer :: EVEC, EVECFF, EVECII
      END MODULE eigv_C
