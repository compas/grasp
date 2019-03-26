!***********************************************************************
!                                                                      *
      SUBROUTINE GETOSD(NAME)
!                                                                      *
!   Interactively determines the data governing the transition prob-   *
!   lem.                                                               *
!                                                                      *
!   Written by Per Jonsson                Last revision: June 1996     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW, NNNP
      USE def_C
      USE decide_C
      USE default_C
      USE foparm_C
      USE grid_C, ONLY: h, hp, n, rnt
      USE npar_C
      USE osc_C, ONLY: LTC
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setiso_I
      USE getrmp_I
      USE setqic_I
      USE radgrd_I
      USE lodrwfi_I
      USE lodrwff_I
      USE brkt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = 128), DIMENSION(2) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      LOGICAL :: YES
      CHARACTER :: CUNITS*4, ANSWER
!-----------------------------------------------
!
!   Determine the physical effects specifications
!
      IF (NDEF /= 0) THEN
         IF (myid .EQ. 0) THEN
            PRINT *, 'The physical speed of light in'
            PRINT *, ' atomic units is',CVAC,';'
            PRINT *, ' revise this value?'
            YES = GETYN ()
            IF (YES) THEN
               PRINT *, 'Enter the revised value:'
               READ *,C
            ELSE
               C = CVAC
            ENDIF
         ENDIF !myid=0
         CALL MPI_Bcast(C,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      ELSE
         C = CVAC
      ENDIF
!
      LFORDR = .FALSE.
      ICCUT = 0
!
!   Determine the multipolarity and parity of the transitions
!
      CALL GETRMP
!
!   Determine the units for the printout and other options
!
      LTC = .FALSE.
!
      IF (NDEF /= 0) THEN
       IF (myid .EQ.0) THEN
         WRITE (6, *) 'Which units are to be used to'
         WRITE (6, *) ' express the transition energies?'
         WRITE (6, *) '    A    : Angstrom:'
         WRITE (6, *) '    eV   : electron volts;'
         WRITE (6, *) '    Hart : Hartree atomic units;'
         WRITE (6, *) '    Hz   : Hertz;'
         WRITE (6, *) '    Kays : Kaysers [cm**(-1)];'
    2    CONTINUE
         READ (*, '(A)') CUNITS
         IF (CUNITS(1:1) == 'A') THEN
            LTC(1) = .TRUE.
         ELSE IF (CUNITS(1:2) == 'eV') THEN
            LTC(2) = .TRUE.
         ELSE IF (CUNITS(1:4) == 'Hart') THEN
            LTC(3) = .TRUE.
         ELSE IF (CUNITS(1:2) == 'Hz') THEN
            LTC(4) = .TRUE.
         ELSE IF (CUNITS(1:4) == 'Kays') THEN
            LTC(5) = .TRUE.
         ELSE
            WRITE (6, *) 'GETOSD: Unable to interpret string;'
            WRITE (6, *) ' reenter ...'
            GO TO 2
         ENDIF
       ENDIF !myid=0
       CALL MPI_Bcast(LTC(1),5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
         LTC(5) = .TRUE.
      ENDIF
!
!      WRITE (6, *) 'Sort transitions by energy?'
!      YES = GETYN()
!      IF (YES) LTC(6) = .TRUE.
!
      IF (NDEF /= 0) THEN
       IF (myid .EQ. 0) THEN
         WRITE (6, *) 'Einstein A and B coefficients are'
         WRITE (6, *) ' printed in SI units; use Hartree'
         WRITE (6, *) ' atomic units instead?'
         YES = GETYN()
         IF (YES) LTC(7) = .TRUE.
       ENDIF !myid=0
       CALL MPI_Bcast(LTC(7),1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
         LTC(7) = .FALSE.
      ENDIF
!
!   Determine the parameters controlling the radial grid
!
      IF (NPARM == 0) THEN
         RNT = EXP((-65.0D00/16.0D00))/Z
         H = 0.5D00**4
         N = MIN(220,NNNP)
      ELSE
!CFF     ... should be Z-dependent
         RNT = 2.0D-06/Z
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D00
      IF (NDEF /= 0) THEN
        IF (myid .EQ. 0) THEN
         WRITE (6, *) 'The default radial grid parameters'
         WRITE (6, *) ' for this case are:'
         WRITE (6, *) ' RNT = ', RNT, ';'
         WRITE (6, *) ' H = ', H, ';'
         WRITE (6, *) ' HP = ', HP, ';'
         WRITE (6, *) ' N = ', N, ';'
         WRITE (6, *) ' revise these values?'
         YES = GETYN()
         IF (YES) THEN
            WRITE (6, *) 'Enter RNT:'
            READ (5, *) RNT
            WRITE (6, *) 'Enter H:'
            READ (5, *) H
            WRITE (6, *) 'Enter HP:'
            READ (5, *) HP
            WRITE (6, *) 'Enter N:'
            READ (5, *) N
         ENDIF
        ENDIF !myid=0
        CALL MPI_Bcast(RNT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast(H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast(HP,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast (N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ENDIF
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!   Set up the coefficients for the numerical proceduRes
!
      CALL SETQIC
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD
!
!   Load the initial state radial wavefunctions.
!
      CALL LODRWFI (NAME(1))
!
!   Load the final state radial wavefunctions.
!
      CALL LODRWFF (NAME(2))

!   Construct the radial overlap matrix.
!
      CALL BRKT

      RETURN
      END SUBROUTINE GETOSD
