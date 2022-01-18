!***********************************************************************
!
      SUBROUTINE SETCSLmpi(NAME, NCORE, IDBLK)
!
!  A container which calls lib92/lodcsh to get
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!  It then calls lodcsl to load the entire file (all blocks), obtaining
!     IQA, JQSA, JCUPA
!  Note that unlike setcsl of mcp[mpi] and rci[mpi]vu, the following
!  quantities are assumed to be known:
!     nblock, ncfblk(), idblk(), ncftot.
!  In this rscfmpi, they have been read into memory from mcp files via
!  a call to setmcp.
!
!  Xinghong He 98-08-06
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE memory_man
      USE hblock_C
      USE def_C
      USE orb_C, NCFTOT=>NCF
      USE stat_C
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodcslmpiGG_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      CHARACTER  :: NAME*(*)
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IOS
      CHARACTER :: STR*15
!-----------------------------------------------
!
! node-0 opens, reads the header of the file

      IF (myid .EQ. 0) THEN
         CALL OPENFL (21, name, 'FORMATTED', 'OLD', ierr)
         IF (ierr .EQ. 1) THEN
            PRINT *, 'Error when opening ',name(1:LEN_TRIM (name))
            STOP
         ENDIF

         READ (21,'(1A15)',IOSTAT = ios) str
         IF ((ios .NE. 0) .OR.                   &
               (str .NE. 'Core subshells:')) THEN
            PRINT *, 'Not a Configuration Symmetry List File;'
            CLOSE (21)
            STOP
         ENDIF

         !..Load header of <name> file
         CALL LODCSH (21, NCORE)
      END IF

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      CALL MPI_Bcast (nw,    1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncore, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nelec, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (np(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nak(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkl(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkj(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nh(1),2*nw, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

! Allocate memories for all blocks and then load the entire file

      CALL ALLOC (iqa, NNNW, NCFTOT, 'IQA', 'SETCSLmpi')
!GG      CALL ALLOC (jqsa, NNNW,3,NCFTOT, 'JQSA', 'SETCSmpiL')
!GG      CALL ALLOC (jcupa, NNNW, NCFTOT, 'JCUPA', 'SETCSLmpi')

!GG      CALL LODCSLmpi (21, NCORE, -119)
      CALL LODCSLmpiGG (21, NCORE, -119)
      ! -119 means "load all blocks"

      RETURN
      END SUBROUTINE SETCSLmpi
