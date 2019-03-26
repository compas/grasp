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
      PROGRAM RSCFVU
!                                                                      *
!   Entry routine for RSCFVU. Controls the entire computation.         *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RSCF92]: CHKPLT, setcsl, setdbg                       *
!                        SETMIX, SETRES, SETSUM, STRSUM.               *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
!   Block version by Xinghong He          Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
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
       USE mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE starttime_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setmcp_I
      USE setcsl_I
      USE getscd_I
      USE strsum_I
      USE setmix_I
      USE factt_I
      USE scf_I
      USE stoptime_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NBLK0 = 50  ! Maximum number of blocks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCORE1, NCOUNT1
      LOGICAL :: EOL, YES
      CHARACTER, DIMENSION(NBLK0) :: IDBLK*8
!-----------------------------------------------
!
!
! Things for timing
!-----------------------------------------------------------------------
      MYID = 0
      NPROCS = 1

      write(*,*)
      write(*,*) 'RMCDHF'
      write(*,*) 'This program determines the radial orbitals   '
      write(*,*) 'and the expansion coefficients of the CSFs         '
      write(*,*) 'in a self-onsistent field proceedure               '
      write(*,*) 'Input file:  isodata, rcsf.inp, rwfn.inp, mcp.30, ...'
      write(*,*)                                                       &
               'Outputfiles: rwfn.out, rmix.out, rmcdhf.sum, rmcdhf.log'
      write(*,*)

      CALL STARTTIME (NCOUNT1, 'RMCDHF')

     OPEN(UNIT=734,FILE='rmcdhf.log',STATUS='UNKNOWN')

!=======================================================================
!  Get NDEF
!=======================================================================

      IF (MYID == 0) THEN
         WRITE (ISTDE, '(A)', ADVANCE='NO') 'Default settings?  (y/n) '
         YES = GETYN()
         IF (YES) THEN
            NDEF = 0
            WRITE(734,'(A)') 'y            ! Default settings'
         ELSE
            NDEF = 1
         ENDIF
      ENDIF

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    SETDBG - Debug output control parameters
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    SETMCP - open and check the  .mcp  files
!    SETCSL - open, check, load data from, and close the  .csl  file
!    STRSUM - append a summary of the inputs to the  .sum  file
!    SETMIX - mixing coefficients file setup
!    FACTT - table of logarithms of factorials setup
!=======================================================================

      CALL SETDBG ('rscf92.dbg')
      CALL SETMC
      CALL SETCON

      CALL SETSUM ('rmcdhf.sum')

      CALL SETMCP (NCORE, NBLK0, IDBLK, 'mcp')
      CALL SETCSL ('rcsf.inp', NCORE1, IDBLK)
      IF (NCORE /= NCORE1) STOP 'rscfvu: ncore'

!=======================================================================
!  Gather all remaining information and perform some setup. This
!  part (routine) asks for user-inputs.
!=======================================================================

      CALL GETSCD (EOL, IDBLK, 'isodata', 'rwfn.inp')

      IF (MYID == 0) THEN
         CALL STRSUM
         IF (EOL) CALL SETMIX ('rmix.out')
      ENDIF

      CALL FACTT

!=======================================================================
!  Proceed with the SCF calculation close all files except
!  the  .sum  file
!=======================================================================

      CALL SCF (EOL, 'rwfn.out')
      CLOSE (734)

!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      CALL STOPTIME (NCOUNT1, 'RMCDHF')

      STOP
      END PROGRAM RSCFVU
