      MODULE setham_I
      INTERFACE
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE SETHAM (myid, nprocs, jblock, ELSTO,ICSTRT, nelmntt,  &
                         atwinv,slf_en)
      USE vast_kind_param,   ONLY: DOUBLE, LONG
      INTEGER(LONG)       :: NELMNTT
      INTEGER             :: JBLOCK, ICSTRT
      INTEGER, INTENT(IN) :: MYID, NPROCS
      REAL(DOUBLE) :: ELSTO, ATWINV
      REAL(DOUBLE), DIMENSION(*) :: slf_en
      END SUBROUTINE SETHAM
      END INTERFACE
      END MODULE
