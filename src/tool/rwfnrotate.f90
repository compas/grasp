!***********************************************************************
!                                                                      *
      PROGRAM RWFNROTATE
!                                                                      *
!   This program rotate orbitals with the same kappa                   *
!                                                                      *
!   Written by Per Jonsson,  Malmo University 1 March 2007             *
!                                                                      *
!***********************************************************************
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
       USE vast_kind_param, ONLY:  DOUBLE
       USE parameter_def
       USE def_C
       USE grid_C
       USE orb_C
       USE wave_C, ONLY: MF, PZ,  PF, QF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setmc_I
      USE setcon_I
      USE setcsla_I
      USE setiso_I
      USE setqic_I
      USE radgrd_I
      USE setrwfa_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER  ::  IFIRST(2), ndef
      LOGICAL  YES
      CHARACTER*24 NAME
      CHARACTER*256 FILNAM
      INTEGER :: i, j, k, MFmax, ncore_not_used
      REAL(double) :: dummy1, dummy2
!

!
      WRITE(*,*)
      WRITE(*,*) 'RWFNROTATE'
      WRITE(*,*) 'This program rotates selected pairs of orbitals'
      WRITE(*,*) 'It is only meaningful to rotate pairs of orbitals'
      WRITE(*,*) 'of the same symmetry'
      WRITE(*,*)
      WRITE(*,*) 'Inputfile : name.w'
      WRITE(*,*) 'Outputfile: name_rot.w'
      WRITE(*,*)

      NDEF = 0

   10 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         PRINT *, 'Names may not start with a blank'
         GOTO 10
      ENDIF
      FILNAM = NAME(1:K-1)//'_rot.w'
      OPEN(36,FILE=FILNAM,FORM='UNFORMATTED',STATUS = 'UNKNOWN')
      WRITE(36) 'G92RWF'

!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Open, check, load data from, and close, the  .csl  file
!
      CALL SETCSLA(NAME,ncore_not_used)
!
!   Read the radial wave functions
!
      CALL GETHFD(NAME)

      WRITE(*,*)
      WRITE(*,*) 'Orbitals'
      WRITE(*,*)
      DO K = 1,NW
         WRITE(*,1000) K, NP(K),NH(K)
      ENDDO

      DO

         WRITE(*,*)
         WRITE(*,*) 'Give the numbers of the two orbitals to be rotated'
         READ(*,*) IFIRST

!
!   Rotate the orbital pair
!

         MFMAX = MAX(MF(IFIRST(1)),MF(IFIRST(2)))
         DO K = 1,MFMAX
            DUMMY1 = PF(K,IFIRST(1))
            DUMMY2 = PF(K,IFIRST(2))
            PF(K,IFIRST(1)) = (DUMMY1 + DUMMY2)/DSQRT(2.d0)
            PF(K,IFIRST(2)) = (DUMMY2 - DUMMY1)/DSQRT(2.d0)

            DUMMY1 = QF(K,IFIRST(1))
            DUMMY2 = QF(K,IFIRST(2))
            QF(K,IFIRST(1)) = (DUMMY1 + DUMMY2)/DSQRT(2.d0)
            QF(K,IFIRST(2)) = (DUMMY2 - DUMMY1)/DSQRT(2.d0)
         ENDDO
         DUMMY1 = PZ(IFIRST(1))
         DUMMY2 = PZ(IFIRST(2))
         PZ(IFIRST(1)) = (DUMMY1 + DUMMY2)/DSQRT(2.d0)
         PZ(IFIRST(2)) = (DUMMY2 - DUMMY1)/DSQRT(2.d0)
         MF(IFIRST(1)) = MFMAX
         MF(IFIRST(2)) = MFMAX

         WRITE(*,*) 'Rotate more orbitals ?'
         YES = GETYN ()
         IF (YES .eqv. .FALSE.) EXIT
      ENDDO


!
      DO I = 1,NW
         WRITE(36) NP(I),NAK(I),E(I),MF(I)
         WRITE(36) PZ(I),(PF(J,I),J = 1,MF(I)),(QF(J,I),J = 1,MF(I))
         WRITE(36) (R(J),J = 1,MF(I))
      ENDDO

      WRITE(*,*)
      WRITE(*,*) 'Execution finished'

1000  FORMAT(I4,I4,A2)
      CONTAINS

!***********************************************************************
!                                                                      *
      SUBROUTINE GETHFD(NAME)
!                                                                      *
!   Interactively determines the data governing the HFS problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE decide_C
      USE def_C
      USE defaulT_C
      USE foparm_C
      USE grid_C,          ONLY: rnt, h, n, hp
      USE npar_C
      USE orb_C
      USE wave_C
      USE wfac_C

      IMPLICIT NONE
!      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
      CHARACTER*24 NAME
!   Open, check, load data from, and close the  .iso  file
!
      CALL SETISO ('isodata')
!
!   Determine the physical effects specifications
!
      IF (NDEF.NE.0) THEN
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
      ELSE
         C = CVAC
      ENDIF
!
!   Determine the parameters controlling the radial grid
!
!   Grid set by SETISO
      IF (NPARM .EQ. 0) THEN
         RNT = EXP (-65.0D00/16.0D00) / Z
         H = 0.5D00**4
         N = MIN (220,NNNP)
      ELSE
         RNT = 2.0D-06
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D00
      IF (NDEF.NE.0) THEN
         PRINT *, 'The default radial grid parameters'
         PRINT *, ' for this case are:'
         PRINT *, ' RNT = ',RNT,';'
         PRINT *, ' H = ',H,';'
         PRINT *, ' HP = ',HP,';'
         PRINT *, ' N = ',N,';'
         PRINT *, ' revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            PRINT *, 'Enter RNT:'
            READ *, RNT
            PRINT *, 'Enter H:'
            READ *, H
            PRINT *, 'Enter HP:'
            READ *, HP
            PRINT *, 'Enter N:'
            READ *, N
         ENDIF
      ENDIF
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD
!
!   Load the radial wavefunctions
!
      CALL SETRWFA(TRIM(NAME)//'.w')
!
      RETURN
      END SUBROUTINE
      END PROGRAM
