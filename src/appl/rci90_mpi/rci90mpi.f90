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
      PROGRAM RCI90mpi
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
      USE mpi_C
      USE decide_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE cslhmpi_I 
      USE setdbg_I 
      USE setmc_I 
      USE setcon_I 
      USE setsum_I 
      USE setres_I 
      USE setmixmpi_I 
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
      INTEGER            :: lenperm, lentmp
      INTEGER :: NCOUNT1, NCOUNT2, NCOUNT_RATE, NCOUNT_MAX 
      INTEGER, DIMENSION(8) :: NYMDUHMSM 
      INTEGER :: K, LENNAME, NCORE, NDUM, J2MAX, NSECONDS 
      LOGICAL :: YES 
      CHARACTER (LEN = 3) :: idstring
!      CHARACTER(LEN=128) :: NAME, TMPDIR, PERMDIR, ISOFILE 
      CHARACTER(LEN=128) :: STARTDIR, NAME, ISOFILE, NAMESAVE,TMPDIR, PERMDIR
      CHARACTER(LEN=8), DIMENSION(NBLK0) :: IDBLK 
      CHARACTER :: CHDATE*8, CHTIME*10, CHZONE*5, STR*8, MSG*128 
!-----------------------------------------------
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================

      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,       &
                           startdir, permdir, tmpdir, 'RCI_MPI')
      WRITE (idstring, '(I3.3)') myid

      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

!=======================================================================
!  Get NDEF on node-0 and then send to all nodes
!=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,*)
         WRITE (ISTDE, *) 'Default settings? ' 
         YES = GETYN() 
         IF (YES) THEN 
            NDEF = 0 
!cjb fort.734
!cjb        write(734,'(A)') 'y            ! Default settings'
         ELSE 
            NDEF = 1 
!cjb fort.734
!cjb        write(734,'(A)') 'n            ! Default settings'
         ENDIF 
      ENDIF 
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     CALL MPI_Bcast (imethod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


!=======================================================================
!  Get name of the state (used in files like <name>.c, <name>.s)
!  info used only on node-0
!=======================================================================

      lenname = 0     ! Otherwise it might cause trouble for myid!=0
      IF (myid .EQ. 0) THEN
         DO
            WRITE (istde,*) 'Name of state: '
            READ (*,'(A)') NAME
            K = INDEX (NAME,' ')
            IF (K .GT. 1) EXIT
            WRITE (istde,*) 'Name may not start with a blank. redo...'
         ENDDO
         NAMESAVE = NAME

         IF (myid .EQ. 0) THEN
            open(unit=734,file=trim(name)//'.clog',status='unknown')
            write(734,'(a)') 'y            ! Default settings'
            write(734,'(a)') trim(name)
         ENDIF

         !...Form the full name of the files used on node-0

         lenname = LEN_TRIM (NAME)
         isofile = permdir(1:lenperm) // '/isodata '
         NAME    = permdir(1:lenperm) // '/' // NAME(1:lenname) // ' '
         lenname = lenperm + lenname + 1
         open(unit=734,file=trim(name)//'.clog',status='unknown')
         write(734,'(a)') 'y            ! Default settings'
         write(734,'(a)') trim(namesave)
      ENDIF

!=======================================================================
   99 CONTINUE


      imcdf = 26        ! Unit for rci.res file
      IPRERUN = 0
!=======================================================================

!=======================================================================
!
!  Checks and settings... Mostly done in backyard.
!
!    SETDBG - Debug output control parameters - all set to false
!    SETMC - machine- and installation- dependent constants
!    SETCON - physical constants
!    SETSUM - open the summary file
!    cslhmpi - load header of the CSL file
!    SETRES - setup the .res  file, one for each node
!    SETMIXmpi - set up mixing coefficients file AND process user input
!    FACTT - table of logarithms of factorials setup
!    SETCSLmpi - open, check, load data from, and close the  .csl  file
!=======================================================================

      CALL SETDBG
      CALL SETMC
      CALL SETCON
      IF (myid .EQ. 0)  CALL SETSUM (NAME)
      CALL cslhmpi (name(1:lenname) // '.c', ncore, nblk0, idblk)
      CALL SETRES (isofile, name(1:lenname) // '.w', idblk)
      CALL SETMIXmpi (NAME, idblk)
      IF (myid .EQ. 0)  CALL STRSUM
      CALL FACTT

!=======================================================================
!   Calculate all the needed Rk integrals - genintrkwrap calls genintrk
!   to compute the integrals (each processor allocates the same amount
!   of memory, but only computes part of the integrals); and then do the
!   "MPI_Allreduce" stuff.
!=======================================================================

      CALL GENINTRKwrap (myid, nprocs, j2max)

!
!   If transverse interaction comput Breit integrals of type 1 and 2
!
      IF (LTRANS) THEN
         CALL GENINTBREIT1WRAP (myid, nprocs, j2max)
         CALL GENINTBREIT2WRAP (myid, nprocs, j2max)
      END IF

!=======================================================================
!   Proceed with the CI calculation
!=======================================================================
      WRITE (6, *) 'Calling MATRIX...' 
 
      CALL MATRIX (NCORE, J2MAX) 
 
      IF (IPRERUN == 1) THEN 
         IPRERUN = 2 
         GO TO 99 
      ENDIF 
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost,          &
                           ncount1, 'RCI_MPI')

      STOP  
      END PROGRAM RCI90mpi
