!***********************************************************************
!                                                                      *
      PROGRAM BIOSCL
!                                                                      *
!     This program calculates the transition parameters for a          *
!     transition between an initial and a final state                  *
!     The program assumes that the radial orbitals of the initial      *
!     and final state have been transformed by the BIOTRA program      *
!     as to become biorthonormal, in which case the normal Racah       *
!     algebra can be used.                                             *
!                                                                      *
!     Written by Per Jonsson,   Department of Computer Science         *
!                             Vanderbilt University, USA               *
!                                                                      *
!     Modified by Gediminas Gaigalas for new spin-angular integration  *
!     and for reducing usage of CPU memory.        NIST, October 2017  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE default_C
      USE debug_C, ONLY: LDBPR, CUTOFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setmc_I
      USE setcon_I
      USE fname_I
      USE mrgcsl_I
      USE setcslm_I
      USE getosd_I
      USE strsum_I
      USE factt_I
      USE oscl_I
      USE mpi_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTEST,NCOUNT1,ILBL,lenperm,lentmp, iii
      LOGICAL :: YES
      CHARACTER, DIMENSION(2) :: NAME*24
      CHARACTER, DIMENSION(2) :: FULLNAME*128
      CHARACTER(LEN=3)  ::  idstring
      CHARACTER(LEN=128) :: ISOFILE,startdir,permdir,tmpdir
!-----------------------------------------------
!
!=======================================================================
!  Start mpi --- get processor info: myid, nprocs, host name; and print
!=======================================================================
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,           &
                      startdir, permdir, tmpdir, 'RTRANSITION_MPI')
      WRITE (idstring, '(I3.3)') myid
      print*, tmpdir , ' = tmpdir'

      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

!=======================================================================
      ISOFILE = trim(startdir)//'/isodata'
      NTEST = 1001

      IF (myid .EQ. 0) THEN
         WRITE (6, *)
         WRITE (6, *) ' Default settings?'
         YES = GETYN()
         WRITE (6, *)
         IF (YES) THEN
            NDEF = 0
            NDUMP = 1
         ELSE
            NDEF = 1
            WRITE (6, *) ' Dump angular data to file?'
            YES = GETYN()
            IF (YES) THEN
               NDUMP = 1
            ELSE
               NDUMP = 0
             ENDIF
         ENDIF
         WRITE (6, *)
         WRITE (6, *) ' Input from a CI calculation?'
         YES = GETYN()
         WRITE (6, *)
         IF (YES) THEN
            INPCI = 0
         ELSE
            INPCI = 1
         ENDIF
      ENDIF   !myid=0
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NDUMP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (INPCI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!Rasa -- start
!GG      LDBPR = .FALSE.
!      WRITE (6, *) ' Generate debug output?'
!      YES = GETYN()
!      WRITE (6, *)
!      IF (YES) THEN
!         LDBPR(18) = .TRUE.
!         WRITE (6, *) ' Enter the cutoff'
!         READ (5, *) CUTOFF
!      ENDIF
!Rasa -- end
!
!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Obtain the names of the initial and final state files
!
      IF (myid .EQ. 0) CALL FNAME(NAME)
      CALL MPI_Bcast (NAME,48,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      FULLNAME(1) = TRIM(startdir) // '/' // TRIM(NAME(1))
      FULLNAME(2) = TRIM(startdir) // '/' // TRIM(NAME(2))
!
!   Open, check, load data from, and close, the initial and final state
!   .csl  files. These files are then merged to one file.
!
!cjb print FULLNAME1= FULLNAME2=
      IF (myid .EQ. 0) then
        print*,"FULLNAME1=",FULLNAME(1)
        print*,"FULLNAME2=",FULLNAME(2)
      end if

      CALL MRGCSL (FULLNAME)
!
!   Open, check, load data from, and close, the merged .csl  file
!
      CALL SETCSLM
!
!   Read mixing coefficients
!
!      CALL READMIX(NAME,INPCI)
!
!   Test mixing coefficients
!
!      IF (NTEST.GT.1000) CALL TESTMIX
!
!   Open, check, load data from, and close the .iso file
!
      CALL SETISO(ISOFILE)

!
!   Get the remaining information
!
      CALL GETOSD (FULLNAME)
!
!   Open and append a summary of the inputs to the  .sum  file
!
      if (myid == 0) then
          iii = len_trim(startdir)
          call sys_chdir(trim(startdir),iii,ierr)
          if (ierr.ne.0) then
              print *, "error changing dir!"
              call exit(1)
          endif
          CALL STRSUM(NAME,INPCI,ILBL)
          iii = len_trim(tmpdir)
          call sys_chdir(trim(tmpdir)//'/'//idstring,iii+4,ierr)
          if (ierr.ne.0) then
              print *, "Error changing dir!"
              call exit(1)
          endif
      endif  !myid = 0
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the transition calculation
!
      CALL OSCL (NAME,FULLNAME,tmpdir,startdir,idstring)
!
!   Print completion message
!
      IF (myid .EQ. 0) THEN
          PRINT *
          PRINT *, 'RTRANSITION_MPI: Execution complete.'
      ENDIF
!=======================================================================
!  Execution finished; Statistics output
!=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost,                     &
                           ncount1, 'RTRANSITION_MPI')
!=======================================================================
!
      STOP
      END PROGRAM BIOSCL
