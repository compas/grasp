!***********************************************************************
!
      SUBROUTINE CSLHMPI(NAME, NCORE, NBLKIN, IDBLK)
!
!  A container which calls setcsll to open, read <name> file to get
!     nblock, ncfblk(), idblk(), ncftot.
!  It then calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  The file pointer points to the first CSL record after this routine.
!
!  Called by rcimpivu, mcpmpi
!
!  Xinghong He 98-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      use mpi_C
      USE HBLOCK_C
      USE def_C, ONLY: EMN,IONCTY,NELEC,Z
      USE ORB_C, NCFTOT=>NCF
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcsll_I
      USE lodcsh_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      INTEGER, INTENT(IN)  :: NBLKIN
      CHARACTER  :: NAME*(*)
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM
!-----------------------------------------------
!

!     node-0 does exactly the same as the serial code does
      IF (myid .EQ. 0) THEN
        CALL ALLOC   (ncfblk, nblkin+1, 'ncfblk', 'cslhmpi')
        CALL SETCSLL (21, name, nblkin, nblock, ncfblk, ncftot, &
                       idblk)
        CALL RALLOC  (ncfblk, nblock+1, 'ncfblk', 'cslhmpi')
         REWIND (21)
         READ   (21,*)
         !..Load header of <name> file
         CALL LODCSH (21, NCORE)
      ENDIF

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      CALL MPI_Bcast (nblock, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncftot, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      IF (myid /= 0) THEN
         CALL alloc (ncfblk, 1+nblock, 'ncfblk', 'cslhmpi')
      ENDIF

      CALL MPI_Bcast (ncfblk, 1+nblock, MPI_INTEGER,0, &
                                              MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (idblk, 8*nblock, MPI_CHARACTER,0, &
                                              MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncore, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nelec, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nw,    1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (np,   nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nak,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkl,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkj,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nh, 2*nw, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE CSLHMPI
