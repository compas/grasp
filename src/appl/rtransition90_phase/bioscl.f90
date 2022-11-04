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
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTEST, NCOUNT1, ILBL
      LOGICAL :: YES
      CHARACTER, DIMENSION(2) :: NAME*24
      CHARACTER(LEN=128) :: ISOFILE
!-----------------------------------------------
!
      ISOFILE = 'isodata'
      NTEST = 1001
      CALL STARTTIME (ncount1, 'RTRANSITION')
      write(*,*)
      write(*,*) 'RTRANSITION'
      write(*,*) 'This program computes transition parameters from'
      write(*,*) 'transformed wave functions'
      write(*,*) 'Input files:  isodata, name1.c, name1.bw, name1.(c)bm'
      write(*,*) '              name2.c, name2.bw, name2.(c)bm         '
      write(*,*) '              optional, name1.lsj.lbl, name2.lsj.lbl'
      write(*,*) '              name1.name2.KT (optional angular files)'
      write(*,*) 'Output files: name1.name2.(c)t                       '
      write(*,*) '              optional, name1.name2.(c)t.lsj         '
      write(*,*) '              name1.name2.KT (angular files)         '
      write(*,*) 'Here K is parity and rank of transition: -1,+1 etc   '




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
!Rasa -- start
      LDBPR = .FALSE.
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
      CALL FNAME (NAME)
!
!   Open, check, load data from, and close, the initial and final state
!   .csl  files. These files are then merged to one file.
!
      CALL MRGCSL (NAME)
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
!   Get the remaining information
!
      CALL GETOSD (ISOFILE,NAME)
!
!   Open and append a summary of the inputs to the  .sum  file
!
      ILBL = 0
      CALL STRSUM (NAME, INPCI,ILBL)
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the transition calculation
!
      CALL OSCL (NAME)
!
!   Print completion message
!
      CALL STOPTIME (ncount1, 'RTRANSITION')
!
      STOP
      END PROGRAM BIOSCL
