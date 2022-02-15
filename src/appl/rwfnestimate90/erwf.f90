
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**               *******  *******   **    **  *******               ***
!**               **       **    **  **    **  **                    ***
!**               **       **    **  **    **  **                    ***
!**               *****    *******   ** ** **  *****                 ***
!**               **       **  **    ** ** **  **                    ***
!**               **       **   **   ** ** **  **                    ***
!**               *******  **    **   **  **   **                    ***
!**                                                                  ***
!**   Program to create a  .rwf  file for RSCF92 by merging one or   ***
!**   more existing  .rwf  files or generating Thomas-Fermi or hy-   ***
!**   drogenic functions.                                            ***
!**                                                                  ***
!**                              GRASP92                             ***
!**           F. A. Parpia, C. F. Fischer, and I. P. Grant           ***
!**                                                                  ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
PROGRAM ERWF
!                                                                      *
!   Entry routine for RCI92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
!                        SETRES, SETSUM, STRSUM.                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
   USE vast_kind_param, ONLY: DOUBLE
   USE DEFAULT_C
   USE CONS_C
   USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   USE getyn_I
   USE setdbg_I
   USE setmc_I
   USE setcon_I
   USE setsum_I
   USE setcsh_I
   USE screenpar_I
   USE getinf_I
   USE strsum_I
   USE genrwf_I
   USE orthsc_I
   USE wrtrwf_I
   IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   INTEGER :: NCORE
   LOGICAL :: YES
!-----------------------------------------------

!   Startup message
!
   WRITE (ISTDE, *) 'RWFNESTIMATE'
   WRITE (ISTDE, *) 'This program estimates radial wave functions'
   WRITE (ISTDE, *) 'for orbitals'
   WRITE (ISTDE, *) 'Input files: isodata, rcsf.inp, optional rwfn file'
   WRITE (ISTDE, *) 'Output file: rwfn.inp'

!
   WRITE (ISTDE, *) 'Default settings ?'
   YES = GETYN()
   IF (YES) THEN
      NDEF = 0
   ELSE
      NDEF = 1
   ENDIF

!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
   CALL SETDBG
!
!   Perform machine- and installation-dependent setup
!
   CALL SETMC
!
!   Set up the physical constants
!
   CALL SETCON
!
!   Open the  .sum  file
!
   IF (NDEF /= 0) CALL SETSUM
!
!   Open, check, load data from, and close the  .csl  file
!
   CALL SETCSH(21, 'rcsf.inp', NCORE)
!
!   Hydrogenic screen parameters for all orbitals
!
   CALL SCREENPAR(NCORE)
!
!   Determine other relevant information
!
   CALL GETINF
!
!   Write the first part of the  .sum  file
!
   IF (NDEF /= 0) CALL STRSUM
!
!   Generate the subshell radial wavefunctions
!
   CALL GENRWF
!
!   Orthogonalize the radial orbitals
!
   CALL ORTHSC
!
!   Write the subshell radial wavefunctions out
!
   CALL WRTRWF
!
!   Print completion message
!
   WRITE (ISTDE, *) 'RWFNESTIMATE: Execution complete.'
!
   STOP
END PROGRAM ERWF
