!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**       ******    *****    *****   *******   *****    *****        ***
!**       **   **  **   **  **   **  **       **   **  **   **       ***
!**       **   **  **       **       **       **   **       **       ***
!**       ******    *****   **       ****      *****       **        ***
!**       **  **        **  **       **          **       **         ***
!**       **   **  **   **  **   **  **         **      **           ***
!**       **   **   *****    *****   **        **      *******       ***
!**                                                                  ***
!**            Relativistic Self-Consistent-Field Program            ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990)                                ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**   Modified by C.F. FIscher for block input                       ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM RSCFmpiVU
!                                                                      *
!   Entry routine for RSCFVU. Controls the entire computation.      *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RSCF92]: CHKPLT, setcsl, setdbg
!                        SETMIX, SETRES, SETSUM, STRSUM.
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
!   Block version by Xinghong He          Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:57:54   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
       USE default_C
       USE core_C
       USE iounit_C
       USE mpi_C
!GG po to isimti
       USE TATB_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbgmpi_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setmcp_I
      USE setcslmpi_I
      USE getscdmpi_I
      USE strsum_I
      USE setmix_I
      USE factt_I
      USE scfmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NBLK0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCORE1, NCOUNT1, lenperm, lentmp
      LOGICAL :: EOL, YES
      CHARACTER, DIMENSION(NBLK0) :: IDBLK*8
!      CHARACTER (LEN = *) :: host
      CHARACTER (LEN = 3) :: idstring
      CHARACTER (LEN = 128) :: startdir,permdir,tmpdir,file_rcsl,file1,file2
!-----------------------------------------------
!
      OPEN(UNIT=734,FILE='rmcdhf.log',STATUS='UNKNOWN')

!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      startdir = '  '  ;    file_rcsl = '  '
      permdir = '  '   ;    file1     = '  '
      tmpdir = '  '    ;    file2     = '  '
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,       &
                      startdir, permdir, tmpdir, 'RMCDHF_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

!=======================================================================
!  Get NDEF
!=======================================================================

      IF (MYID == 0) THEN
         WRITE (istde,*)
         WRITE (ISTDE, '(A)') 'Default settings?  (y/n) '
         YES = GETYN()
         IF (YES) THEN
            NDEF = 0
!cjb fort.734
!cjb        WRITE(734,'(A)') 'y            ! Default settings'
         ELSE
            NDEF = 1
         ENDIF
      ENDIF
     CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    CHKPLT - compatibility of plant substitutions
!    SETDBGmpi - Debug output control parameters
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    SETMCP - open and check the  .mcp  files
!    SETCSLmpi - open, check, load data from, and close the  .csl  file
!    STRSUM - append a summary of the inputs to the  .sum  file
!    SETMIX - mixing coefficients file setup
!    FACTT - table of logarithms of factorials setup
!=======================================================================
      CALL SETDBGmpi (permdir(1:lenperm) // '/rscf92.dbg')
      CALL SETMC
      CALL SETCON

      IF (myid .EQ. 0) CALL SETSUM (permdir(1:lenperm) // '/rmcdhf.sum')

      CALL SETMCP (NCORE, NBLK0, IDBLK, 'mcp' // idstring)
      if(myid == 0) file_rcsl = permdir(1:lenperm)//'/rcsf.inp'
      CALL SETCSLmpi (file_rcsl, ncore1, idblk)
      IF (NCORE /= NCORE1) STOP 'rscfmpivu: ncore'

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      if(myid == 0) file1 = permdir(1:lenperm) // '/isodata'
      if(myid == 0) file2 = permdir(1:lenperm) // '/rwfn.inp'
      CALL GETSCDmpi (EOL, idblk, file1, file2 )
!     &            permdir(1:lenperm) // '/isodata',
!     &            permdir(1:lenperm) // '/rwfn.inp')
      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm) // '/rmix.out'
      IF (myid .EQ. 0) THEN
         CALL STRSUM
         IF (EOL) CALL SETMIX (file1)
      ENDIF

      CALL FACTT

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm) // '/rwfn.out'

      CALL scfmpi (EOL, file1)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost, ncount1, 'RMCDHF_MPI')

      END PROGRAM RSCFmpiVU
