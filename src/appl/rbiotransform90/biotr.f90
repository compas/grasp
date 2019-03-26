!***********************************************************************
!                                                                      *
!     This program transforms the initial and final state radial       *
!     orbitals to a biorthonormal set. The CI coefficients are         *
!     then counter transformed as to leave both the initial and        *
!     final state wave functions invariant                             *
!                                                                      *
!     General references:                                              *
!                                                                      *
!     P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)  *
!     J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and         *
!     C. Froese Fischer, Phys. Rev. E, 4499 (1995)                     *
!                                                                      *
      PROGRAM BIOTR
!                                                                      *
!     Program written by                                               *
!                                                                      *
!     Per Jonsson, Department of Computer Science                      *
!                  Vanderbilt University, USA                          *
!                                                                      *
!     Date June 1996                                                   *
!                                                                      *
!     Modified by Gediminas Gaigalas for new spin-angular integration  *
!     and for reducing usage of CPU memory.        NIST, October 2017  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:17:22   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW, NNNP
      USE default_C,       ONLY: ndef, ndump
      USE sbdat_C,         ONLY: NLMAX, KAMAX, NSHLII, NSHLFF
      USE biorb_C,         ONLY: nwii, nwff
      USE wave_C,          ONLY: pfii, qfii, pfff, qfff
      USE cimat_C,         ONLY: cici, cfci
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setmc_I
      USE setcon_I
      USE setiso_I
      USE radpar_I
      USE radgrd_I
      USE setqic_I
      USE fname_I
      USE setcslb_I
      USE tcsl_I
      USE kapdata_I
      USE lodrwfi_I
      USE lodrwff_I
      USE brkt_I
      USE gets_I
      USE biotr1_I
      USE radfile_I
      USE genmcp_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!GG      INTEGER, PARAMETER :: LWORK1 = 100000
      INTEGER, PARAMETER :: LWORK1 = 10000000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NLMAX) :: NINSHLI, NINSHLF
      INTEGER :: NTESTG, NTESTL, NTEST, INPCI, NCORE1, NCORE2, MXL, K
      INTEGER :: NCOUNT1
      REAL(DOUBLE), DIMENSION(LWORK1) :: WORK
      REAL(DOUBLE), DIMENSION(NLMAX*NLMAX) :: CISHL
      REAL(DOUBLE), DIMENSION(NNNW,NNNW)   :: S
      REAL(DOUBLE), DIMENSION(NLMAX*NLMAX) :: CFSHL
      LOGICAL :: YES
      CHARACTER, DIMENSION(2) :: NAME*24
      CHARACTER(LEN=128) :: ISOFILE
!-----------------------------------------------
!
      ISOFILE = 'isodata'
!   Debug flags
!
      NTESTG = 0
      NTESTL = 0
      NTEST = MAX0(NTESTL,NTESTG)
!
      CALL STARTTIME (ncount1, 'BIOTRANSFORM')
      write(*,*)
      write(*,*) 'RBIOTRANSFORM'
      write(*,*) 'This program transforms the initial and final wave'
      write(*,*) 'functions so that standard tensor albegra can be'
      write(*,*) 'used in evaluation of the transition parameters'
      write(*,*) 'Input files:  isodata, name1.c, name1.w, name1.(c)m'
      write(*,*) '              name2.c, name2.w, name2.(c)m         '
      print*,'              name1.TB, name2.TB (optional angular files)'
      write(*,*) 'Output files: name1.bw, name1.(c)bm, '
      WRITE(*,*) '              name2.bw, name2.c(bm)'
      write(*,*) '              name1.TB, name2.TB (angular files)'
      write(*,*)

      WRITE (6, *) 'Default settings?'
      YES = GETYN()
      WRITE (6, *)
      IF (YES) THEN
         NDEF = 0
         NDUMP = 1
      ELSE
         NDEF = 1
         WRITE (6, *) 'Dump angular data on file?'
         YES = GETYN()
         WRITE (6, *)
         IF (YES) THEN
            NDUMP = 1
         ELSE
            NDUMP = 0
         ENDIF
      ENDIF
      WRITE (6, *) 'Input from a CI calculation?'
      YES = GETYN()
      WRITE (6, *)
      IF (YES) THEN
         INPCI = 0
      ELSE
         INPCI = 1
      ENDIF
!
!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Open, check, load data from, and close the  .iso  file
!
      CALL SETISO (ISOFILE)
!
!   Determine the parameters controlling the radial grid
!
      CALL RADPAR
!
!   Generate the radial grid
!
      CALL RADGRD
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Obtain the names of the initial and final state files
!   and open files where the transformed orbitals and CI
!   coefficients are to be dumped
!
      CALL FNAME (NAME)
!
!   Open, check, load data from and close the initial state CSL file.
!
      CALL SETCSLB (NAME(1), NCORE1,1)
!
!   Transfer the data to the initial state COMMON
!
      CALL TCSL (1)
!
!   Open, check, load data from and close the final state CSL file.
!
      CALL SETCSLB (NAME(2), NCORE2,2)
!
!   Transfer the data to the final state COMMON
!
      CALL TCSL (2)
!
!   Determine the number of kappa quantum numbers and
!   the number of orbitals for each kappa quantum number
!   for the initial state and final states
!
      CALL KAPDATA (NTESTG, NCORE1, NCORE2)
!
!   Read the the radial orbitals for the initial state
!
      CALL LODRWFI (NAME(1), NTESTG)
!
!   Read the the radial orbitals for the initial state
!
      CALL LODRWFF (NAME(2), NTESTG)
!
!   Calculate the radial overlap matrices
!
      WRITE (*, *)
      WRITE (*, *) ' ******************************************'
      WRITE (*, *) '  Overlap matrix before orbital rotations'
      WRITE (*, *) ' *****************************************'
      WRITE (*, *)

      CALL BRKT

      CALL GETS (S, NWII, NWFF)

!
!   Once we have the overlap matrices
!   we can manipulate the initial and final state separately.
!
      MXL = KAMAX
!
!. Calculate biorthonormal orbitals, and orbital matrix
!. for counter transformation of CI coefficients.
!
      CALL BIOTR1 (PFII, QFII, NSHLII, NINSHLI, PFFF, QFFF, &
         NSHLFF, NINSHLF, NNNP, KAMAX, WORK, LWORK1, NTESTG, &
         CISHL, CICI, CFSHL, CFCI)
      WRITE (*, *)
      WRITE (*, *) ' ****************************************'
      WRITE (*, *) '  Overlap matrix after orbital rotations'
      WRITE (*, *) ' ****************************************'
      WRITE (*, *)

      CALL BRKT
!
!  Write the transformed radial functions to file
!
      CALL RADFILE (NAME)
!
!   Obtain one-electron coupling coefficients for the initial state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      CALL GENMCP (NAME(1), 1, NTESTG, INPCI)
!***** added by Yu Zou, Feb.18,2000 ***********
      DO K = 1, KAMAX
         CLOSE(UNIT=80 + K, STATUS='DELETE')
      END DO
!***** added by Yu Zou, Feb.18,2000 ***********
!
!   Obtain one-electron coupling coefficients for the final state.
!   The coefficients are dumped on files one kappa in turn and
!   thus the different kappa can be manipulated independently.
!   The interface with the transformation part is in the routine mcp
!
      CALL GENMCP (NAME(2), 2, NTESTG, INPCI)

      DO K = 1, KAMAX
         CLOSE(UNIT=80 + K, STATUS='DELETE')
      END DO
      CALL STOPTIME (ncount1, 'BIOTRANSFORM')

      STOP
      END PROGRAM BIOTR
