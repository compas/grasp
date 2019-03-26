
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
      PROGRAM GENMCP
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
!   Updated by Charlotte F. Fischer                                    *
!                                                                      *
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE default_C
      Use hblock_C
      Use orb_C
      Use stat_C
      Use mcp_C
      Use iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE starttime_I
      USE setdbg_I
      USE setmc_I
      USE setsum_I
      USE cslh_I
      USE setmcp2_I
      USE strsum_I
      USE factt_I
!      USE settmp_I
      USE settmpgg_I
      USE stoptime_I
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
      CHARACTER(LEN=128) :: STARTDIR*128, PERMDIR*128, TMPDIR*128
      CHARACTER(LEN=8), DIMENSION(NBLK0):: IDBLK
      INTEGER :: MYID,NPROCS
!-----------------------------------------------------------------------
!
      write(*,*)
      write(*,*) 'RANGULAR'
      write(*,*) 'This program performs angular integration '
      write(*,*) 'Input file:  rcsf.inp'
      write(*,*) 'Outputfiles: mcp.30, mcp.31, ....'
      write(*,*) 'rangular.log                     '
      write(*,*)
      OPEN(UNIT=739,FILE='rangular.log',STATUS='UNKNOWN')
      MYID = 0
      NPROCS = 1
      CALL STARTTIME (NCOUNT1, 'RANGULAR')
!=======================================================================
!  Get NDEF
!=======================================================================
      IF (MYID == 0) THEN
!         WRITE (ISTDE, '(A)', ADVANCE='NO') 'Default settings?  (y/n) '
         WRITE (istde,'(A)') ' Full interaction?  (y/n) '
         YES = GETYN()
         IF (YES) THEN
            NDEF = 0
            write(739,'(A)') 'y            ! Full interaction'
         ELSE
            NDEF = 1
            write(739,'(A)') 'n            ! Full interaction'
         ENDIF
      ENDIF
!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    chkplt - compatibility of plant substitutions
!    setdbg - debug output control parameters
!    setmc - machine- and installation- dependent constants
!    setsum - open the summary file
!    cslh - load header of the csl file
!    setmcp - open and check the  .mcp  files
!    strsum - append a summary of the inputs to the  .sum  file
!    factt - table of logarithms of factorials setup
!=======================================================================
!CFF  delete chkplt
!     CALL CHKPLT ('GENMCP')
      CALL SETDBG (DEBUG, 'genmcp.dbg')
      CALL SETMC
!      IF (NDEF/=0 .AND. MYID==0) CALL SETSUM ('genmcp.sum')
      CALL CSLH ('rcsf.inp', NCORE, NBLK0, IDBLK)
      RESTRT = .FALSE.
      CALL SETMCP2 (MYID, NPROCS, NCORE, IDBLK, 'mcp')
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
         CALL LODCSH2 (21, NCORE, NB)
         !*** Open tmp.xx files for block nb ***
!GG         CALL SETTMP (NB, KMAX, 'tmp')
         CALL SETTMPGG (nb, 30, 'tmp')
         !*** Generation of MCP coefficients ***
         CALL MCP (NB, RESTRT, MYID, NPROCS, 'mcp')
      END DO
      CLOSE(24)                                  ! Summary file
      CLOSE(739)                                  ! rangular.log
      IF (DEBUG) CLOSE(99)                       ! Debug file
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================
      CALL STOPTIME (NCOUNT1, 'RANGULAR')
      STOP
      END PROGRAM GENMCP
