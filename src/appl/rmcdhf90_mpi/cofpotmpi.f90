      SUBROUTINE COFPOTmpi(EOL, J, NPTS)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE pote_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcof_I
      USE ypot_I
      USE xpot_I
      USE lagcon_I
      USE dacon_I
      USE POTE_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      INTEGER  :: NPTS
      LOGICAL  :: EOL
      REAL(DOUBLE), DIMENSION(NPTS) :: TMPMPI
!-----------------------------------------------------------------------

      CALL SETCOF (EOL, J)
      CALL YPOT (J)
      CALL XPOT (J)
      CALL LAGCON (J, NPROCS)
      CALL DACON

      CALL MPI_Allreduce (YP, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, YP, 1)

      CALL MPI_Allreduce (XP, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, XP, 1)

      CALL MPI_Allreduce (XQ, tmpmpi, npts, MPI_DOUBLE_PRECISION,  &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (npts, tmpmpi, 1, XQ, 1)

      RETURN
      END SUBROUTINE COFPOTmpi
