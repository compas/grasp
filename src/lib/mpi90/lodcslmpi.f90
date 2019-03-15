!***********************************************************************
      SUBROUTINE lodcslmpi (nfile, ncore, jblock)

! An MPI container of lodcsh2 which loads CSL list of the current block
! into memory. It forwards the call together with the same set of
! parameters to lodcsh2 and then broadcasts the results to all nodes.
!
! Note: Memories have been allocated/deallocated each block outside.
! This subroutine calls lodcsh2 on node-0 to generate the data for the
! block; and then broadcasts to all other nodes. A new MPI data type
! of 4 byte-long is created to handle 64-bit machines whose MPI
! implementation does not support 4-byte integers. If jblock=-119,
! then ALL blocks will be loaded instead of just one. This is
! implemented in lodcsh2.
!
! Currently used by rcimpivu, mcpmpi, rscfmpivu
!
! Xinghong He 98-08-06
!
!***********************************************************************
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  BYTE
      USE parameter_def,   ONLY:  NNNW
      use mpi_C
      USE ORB_C,           ONLY: NCF, IQA
      USE STAT_C,          ONLY: JQSA, JCUPA

      IMPLICIT NONE
      INTEGER       :: nfile, ncore, jblock
      INTEGER(BYTE) :: MPIX_INT1

!-----------------------------------------------------------------------

      IF (myid .EQ. 0) THEN
         CALL lodcsh2 ((nfile), (ncore), (jblock))
      ENDIF

! Construct mpi data type for Integer*1 and then broadcast.

!cjb mpix_bytes
!cjb if your compiler or your MPI does not accept type MPIX_INT1
!cjb try standard MPI type MPI_INTEGER1
      CALL mpix_bytes (1, MPIX_INT1, ierr)
!     CALL MPI_Type_create_f90_integer (15, MPIX_INT1, ierr)

      CALL MPI_Bcast                                                   &
               (IQA(:,:),     NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!cjb mpix_bytes
!cjb           (IQA(:,:),    NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast                                                   &
               (JQSA(:,:,:),3*NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!cjb mpix_bytes
!cjb           (JQSA(:,:,:),3*NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast                                                   &
               (JCUPA(:,:),   NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!cjb mpix_bytes
!cjb           (JCUPA(:,:),   NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE lodcslmpi
