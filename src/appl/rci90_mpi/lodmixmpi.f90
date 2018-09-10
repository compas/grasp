!***********************************************************************
!                                                                      *
      SUBROUTINE LODMIXmpi (idblk)
!                                                                      *
!   Determines the eigenpairs required;  this information is written   *
!   to the head of the  .mix  file.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, lodstate                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Nov 1992   *
!   Block version by Xinghong He          Last revision:  9 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE def_C
      USE hblock_C
      USE blk_C,           ONLY: idxblk
      USE orb_C,           ONLY: ncf, nw, iqa
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lodstate_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=8),DIMENSION(*), INTENT(IN):: idblk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ncftot, noffset, nvecsiz, jb, j
      INTEGER :: NELEC4, ncftot4, NW4, ncmin4, nvecsiz4, nblock4
!-----------------------------------------------
 
! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here
 
      CALL alloc (ncmaxblk, nblock, 'NCMAXBLK', 'LODMIXmpi')
      CALL alloc (nevblk, nblock, 'NEVBLK', 'LODMIXmpi')
 
!      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
!      stop
      IF (myid .EQ. 0) THEN
         CALL LODSTATE (IDBLK)
      ENDIF

      CALL MPI_Bcast (nevblk, nblock, MPI_INTEGER, 0,               &
                                       MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ncmaxblk, nblock, MPI_INTEGER, 0,             &
                                       MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ncmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      IF (ncmin .NE. 0) THEN
         IF (myid .NE. 0) THEN   ! lodstate allocated it from node-0
            CALL alloc (iccmin, ncmin, 'ICCMIN', 'LODMIXmpi')
         ENDIF

         CALL MPI_Bcast (iccmin, ncmin, MPI_INTEGER4, 0,            &
                                       MPI_COMM_WORLD, ierr)
      ENDIF

!
!   Determine other auxiliary arrays, parameters
!
!    idxblk() is the block number of an eigenstate
!    nvecsiz  is the total size of the eigenvector array
!
      CALL alloc (idxblk, ncmin, 'IDXBLK', 'LODMIXmpi')
      ncftot = 0
      noffset = 0
      nvecsiz = 0
      DO jb = 1, nblock
         DO j = 1, nevblk(jb)
            idxblk(j + noffset) = jb
         ENDDO
         ncftot = ncftot + ncfblk(jb)
         noffset = noffset + nevblk(jb)
         nvecsiz = nvecsiz + ncfblk(jb) * nevblk(jb)
      ENDDO
      IF (noffset .NE. ncmin) STOP 'lodmixmpi: ncmin trouble'
!
!   Write header data to the  .mix  file
!
      IF (myid .EQ. 0) THEN
         NELEC4 = NELEC
         ncftot4 = ncftot
         NW4 = NW
         ncmin4 = ncmin
         nvecsiz4 = nvecsiz
         nblock4 = nblock
         WRITE (25) NELEC4, ncftot4, NW4, ncmin4, nvecsiz4, nblock4
      ENDIF
!
      RETURN
      END SUBROUTINE LODMIXmpi
