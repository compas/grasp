      MODULE iniestsd_I
      INTERFACE
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE INIESTSD (nmax,ncf, myid,nprocs,NIV,BASIS,IMCDF,EAV)
      USE vast_kind_param,ONLY: DOUBLE
      INTEGER, INTENT(IN) :: nmax, ncf, myid, nprocs, niv, imcdf
      REAL(DOUBLE), INTENT(IN) :: EAV
      REAL(DOUBLE), DIMENSION(*), INTENT(IN):: Basis
      END SUBROUTINE
      END INTERFACE
      END MODULE
