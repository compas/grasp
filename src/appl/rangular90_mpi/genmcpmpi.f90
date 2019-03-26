!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**        *****   ******  **   **  **    **   *****   ******        ***
!**       **   **  **      ***  **  ***  ***  **   **  **   **       ***
!**       **       **      ***  **  ** ** **  **       **   **       ***
!**       **  ***  ****    ** ****  ** ** **  **       ******        ***
!**       **   **  **      **  ***  **    **  **       **            ***
!**       **   **  **      **   **  **    **  **   **  **            ***
!**        *****   ******  **   **  **    **   *****   **            ***
!**                                                                  ***
!**   Program for generating the energy expression for H(DC). This   ***
!**   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
!**   and C. F. Fischer, 1990).                                      ***
!**                                                                  ***
!**                         GRASP92 Version                          ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM GENMCPMPI
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer
!                                                                      *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mpi_C
      USE default_C
      Use hblock_c
      Use orb_C
      Use stat_c
      Use mcp_c
      Use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setsum_I
      USE cslh_I
      USE strsum_I
      USE factt_I
      USE settmp_I
      USE lodcslmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER NBLK0, NCOUNT1, NCORE, NB
      PARAMETER (NBLK0 = 50)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  DEBUG, RESTRT, YES
      CHARACTER(LEN=128) :: STARTDIR, PERMDIR, TMPDIR, FILE_RCSL
      CHARACTER(LEN=8), DIMENSION(NBLK0):: IDBLK
      CHARACTER(LEN=3) :: IDSTRING
      INTEGER :: LENPERM, LENTMP
!-----------------------------------------------
      OPEN(UNIT=739,FILE='rangular.log',STATUS='UNKNOWN')
      startdir = ' '
      permdir = ' '
      tmpdir = ' '
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, 'RANGULAR_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)
!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================
      IF (myid == 0) THEN
         WRITE (istde,*)
!         WRITE (istde,'(A)') ' Default settings?  (y/n) '
         WRITE (istde,'(A)') ' Full interaction?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
            write(739,'(A)') 'y            ! Full interaction'
         ELSE
            NDEF = 1
            write(739,'(A)') 'n            ! Full interaction'
         ENDIF
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================

      CALL setdbgmpi (DEBUG, permdir(1:lenperm) // '/genmcp.dbg')
      CALL SETMC
      if(myid == 0) file_rcsl = permdir(1:lenperm)//'/rcsf.inp'
      CALL cslhmpi (file_rcsl, ncore, nblk0, idblk)
      RESTRT = .FALSE.
!cjb  myid, nprocs = NOT args
!cjb  CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)
      CALL setmcpmpi (ncore, idblk, 'mcp' // idstring)
      IF (NDEF/=0 .AND. MYID==0) CALL STRSUM
      CALL FACTT

!=======================================================================
!     For each block, generate and sort the data
!=======================================================================

      DO NB = 1, NBLOCK
         NCF = NCFBLK(NB)                        ! This ncf goes to common
         IF (MYID == 0) THEN
            WRITE (6, *)
            WRITE (6, *) 'Block ', NB, ',  ncf = ', NCF
         ENDIF
            !*** Load current CSL block. Memories de-allocated in mcp ***
         CALL ALLOC (iqa, NNNW, NCF, 'IQA', 'GENMCP')
         CALL ALLOC (jqsa, NNNW,3,NCF, 'JQSA', 'GENMCP')
         CALL ALLOC (jcupa, NNNW, NCF, 'JCUPA', 'GENMCP')
!
         CALL LODCSLmpi (21, NCORE, NB)
         !*** Open tmp.xx files for block nb ***
         CALL SETTMP (NB, KMAX, 'tmp' // idstring)
         !*** Generation of MCP coefficients ***
!cjb  myid, nprocs = NOT args
!cjb     CALL MCPmpi (NB, RESTRT, MYID, NPROCS, 'mcp' // idstring)
         CALL MCPmpi (NB, RESTRT, 'mcp' // idstring)
      END DO
      CLOSE(24)                                  ! Summary file
      CLOSE(739)                                 ! rangular.log
      IF (DEBUG) CLOSE(99)                       ! Debug file
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      CALL stopmpi2 (myid,nprocs,host,lenhost,ncount1,'RANGULAR_MPI')
      STOP
      END PROGRAM GENMCPMPI
