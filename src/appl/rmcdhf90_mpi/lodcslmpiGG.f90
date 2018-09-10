!***********************************************************************
      SUBROUTINE lodcslmpiGG (nfile, ncore, jblock)

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
      USE vast_kind_param, ONLY: BYTE
      USE parameter_def,   ONLY: NNNW
      USE ORB_C,           ONLY: NCF, IQA
      USE syma_C,          ONLY: JPGG, nblk0
      use mpi_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: nfile, ncore, jblock
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER(BYTE) :: MPIX_INT1
!-----------------------------------------------------------------------
!
      IF (myid == 0) THEN
!GG         CALL lodcsh2 ((nfile), (ncore), (jblock))
         CALL lodcsh2GG ((nfile), (ncore), (jblock))
      ENDIF

! Construct mpi data type for Integer*1 and then broadcast.
!cjb mpix_bytes
!cjb if your compiler or your MPI does not accept type MPIX_INT1
!cjb try standard MPI type MPI_INTEGER1
      CALL mpix_bytes (1, MPIX_INT1, ierr)
      CALL MPI_Bcast (IQA(:,:),NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!cjb  CALL MPI_Bcast (IQA(:,:),NNNW*NCF,MPI_INTEGER1,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (JPGG(:),nblk0,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!GG      CALL MPI_Bcast                                                   &
!GG               (JQSA(:,:,:),3*NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
!GG      CALL MPI_Bcast                                                   &
!GG               (JCUPA(:,:),   NNNW*NCF,MPIX_INT1,0,MPI_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE lodcslmpiGG
