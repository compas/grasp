!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**                                                                  ***
!**             ******    *****   ****   *****    *****              ***
!**             **   **  **   **   **   **   **  **   **             ***
!**             **   **  **        **   **   **       **             ***
!**             ******   **        **    *****       **              ***
!**             **  **   **        **      **       **               ***
!**             **   **  **   **   **     **      **                 ***
!**             **   **   *****   ****   **      *******             ***
!**                                                                  ***
!**          Relativistic Configuration-Interaction Program          ***
!**                                                                  ***
!**   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
!**   Grant, and C. F. Fischer, 1990).                               ***
!**                                                                  ***
!**                            GRASP92                               ***
!**          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
      PROGRAM RCI92 
!                                                                      *
!   Entry routine for RCI92. Controls the entire computation.          *
!                                                                      *
!   Call(s) to: [LIB92]: SETMC, SETCON.                                *
!               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
!                        SETRES, SETSUM, STRSUM.                       *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Updated by Xinghong He                Last revision: 23 Jun 1998   *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE default_C
      USE blim_C
      USE where_C
      USE cons_C
      USE def_C
      USE hblock_C
      USE blk_C
      USE iccu_C
      USE iounit_C
      USE terms_C
      USE decide_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE setdbg_I 
      USE setmc_I 
      USE setcon_I 
      USE setsum_I 
      USE setcsl_I 
      USE setres_I 
      USE setmix_I 
      USE strsum_I 
      USE factt_I 
      USE matrix_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NBLK0 = 50 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: nblk=50
      INTEGER :: NCOUNT1, NCOUNT2, NCOUNT_RATE, NCOUNT_MAX 
      INTEGER, DIMENSION(8) :: NYMDUHMSM 
      INTEGER :: MYID, NPROCS, K, LENNAME, NCORE, NDUM, J2MAX, NSECONDS 
      LOGICAL :: YES 
      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE 
      CHARACTER(LEN=8), DIMENSION(NBLK0) :: IDBLK 
      CHARACTER :: CHDATE*8, CHTIME*10, CHZONE*5, STR*8, MSG*128 
!-----------------------------------------------
!
 
      IMCDF = 26                                 ! Unit for rci.res file 
      IPRERUN = 0 
      MYID = 0 
      NPROCS = 1 

      write(*,*)
      write(*,*) 'RCI'
      write(*,*) 'This is the configuration interaction program '
      write(*,*) 'Input file:  isodata, name.c, name.w'
      write(*,*) 'Outputfiles: name.cm, name.csum, name.clog '
      write(*,*) '             rci.res (can be used for restart)'
      write(*,*)

 
!
! Start timing
!
      CALL STARTTIME (ncount1, 'RCI')
!
! Get NDEF
!
!      WRITE (ISTDE, *) 'RCI2: Execution begins ...' 
!      WRITE (ISTDE, *) 
      WRITE (ISTDE, '(A)') 'Default settings? ' 
      YES = GETYN() 
      IF (YES) THEN 
         NDEF = 0 
      ELSE 
         NDEF = 1 
      ENDIF 
!
! Get name of the state (used in files like <name>.c, <name>.s)
!
      DO WHILE(.TRUE.) 
         WRITE (ISTDE, '(A)') 'Name of state: ' 
         READ (*, '(A)') NAME 
         K = INDEX(NAME,' ') 
         IF (K > 1) EXIT  
         WRITE (ISTDE, *) 'Name may not start with a blank. redo...' 
      END DO 

! Now the name of the state is known, open the log file

      open(unit=734, file=trim(name)//'.clog',status='unknown')
      write(734,'(a)') 'y            ! Default settings'
      write(734,'(a)') trim(name)
 
!         ...Form the full name of the files used on node-0
 
      LENNAME = LEN_TRIM(NAME) 
      ISOFILE = 'isodata' 
!      WRITE (6, *) 'isofile = ', ISOFILE(1:LEN_TRIM(ISOFILE)) 
!      WRITE (6, *) 'name = ', NAME(1:LEN_TRIM(NAME)) 

   99 CONTINUE 
 
!
! In SETDBG of this version all control logicals are set to
! false thus no debug output will be made
!
!      WRITE (6, *) 'Calling SETDBG...' 
      CALL SETDBG 
!
! Perform machine- and installation-dependent setup
!
!      WRITE (6, *) 'Calling SETMC...' 
      CALL SETMC 
!
! Set up the physical constants
!
!      WRITE (6, *) 'Calling SETCON...' 
      CALL SETCON 
!
! Open summary file
!
!      WRITE (6, *) 'Calling SETSUM...' 
      CALL SETSUM (NAME) 
 
!      WRITE (6, *) 'Calling setcsl...' 
      CALL SETCSL (NAME(1:LENNAME)//'.c', NCORE, NBLK0, IDBLK) 
!
! Set up the  .res  file; determine if this is a restart.
!
!      WRITE (6, *) 'Calling SETRES...' 
      CALL SETRES (ISOFILE, NAME(1:LENNAME)//'.w', IDBLK) 
!
!   Open the  .mix  file; determine the eigenpairs required
!
!      WRITE (6, *) 'Calling SETMIX...' 
      CALL SETMIX (NAME, IDBLK) 
!
!   Append a summary of the inputs to the  .sum  file
!
      WRITE (6, *) 'Calling STRSUM...' 
      CALL STRSUM 
!
!   Set up the table of logarithms of factorials
!
      WRITE (6, *) 'Calling FACTT...' 
      CALL FACTT 
!
!   Calculate all the needed Rk integrals
!
      WRITE (6, *) 'Calling GENINTRK...' 
      CALL GENINTRK (MYID, NPROCS, NDUM, J2MAX) 
!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      IF (LTRANS) THEN
         PRINT *, 'Calling GENINTBREIT1...'
         CALL GENINTBREIT1 (myid, nprocs, ndum, j2max)
         PRINT *, 'Calling GENINTBREIT2...'
         CALL GENINTBREIT2 (myid, nprocs, ndum, j2max)
      END IF
!
!   Proceed with the CI calculation
!
      WRITE (6, *) 'Calling MATRIX...' 
 
      CALL MATRIX (NCORE, J2MAX) 
 
      IF (IPRERUN == 1) THEN 
         IPRERUN = 2 
         GO TO 99 
      ENDIF 
 
      IF (MYID == 0) THEN 
         WRITE (6, *) 
         WRITE (6, *) 
         WRITE (6, *) 'Finish time, Statistics' 
         WRITE (6, *) 
      ENDIF 
      close(734)
      CALL STOPTIME (ncount1, 'RCI') 
!      CALL SYSTEM_CLOCK (NCOUNT2, NCOUNT_RATE, NCOUNT_MAX) 
!      NCOUNT2 = NCOUNT2 - NCOUNT1 
!      NSECONDS = NCOUNT2/NCOUNT_RATE 
!      WRITE (STR, '(I8)') NSECONDS 
!      MSG = STR//' seconds ' 
!      WRITE (6, *) MSG 
! 
!      CALL DATE_AND_TIME (CHDATE, CHTIME, CHZONE, NYMDUHMSM) 
! 
!      MSG = ' Date: '//CHDATE//' Time: '//CHTIME//' Zone: '//CHZONE 
!      WRITE (6, *) MSG 
!
!   Print completion message
!
!      WRITE (6, *) 'RCI2: Execution complete.' 
!
      STOP  
      END PROGRAM RCI92 
