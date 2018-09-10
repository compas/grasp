!***********************************************************************
!                                                                      *
      SUBROUTINE GETCID (isofile, rwffile, idblk)
!                                                                      *
!   Interactively determines the data governing the CI problem.        *
!   iccut is replaced by an array iccutblk(1:nblock)
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE decide_C
      USE def_C,           ONLY: cvac, c, z, accy, nelec, emn
      USE default_C
      USE grid_C
      USE hblock_C
      USE iccu_C
      USE iounit_C
      USE npar_C
      USE npot_C,          ONLY: zz, nnuc
      USE nsmdat_C
      USE orb_C
      USE wave_C
      USE wfac_C
      USE blim_C
      USE where_C
      USE qedcut_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setisompi_I
      USE setrwfmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      CHARACTER(LEN=*) :: isofile, rwffile
      CHARACTER(LEN=8), DIMENSION(*) ::  idblk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL :: getyn, yes
       INTEGER :: i, jblock, ntmp, np10, j
!-----------------------------------------------------------------------
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
! common block.
!
      PRINT *, 'Calling SETISO ...'
      CALL SETISOmpi (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      IF (NDEF .NE. 0) THEN
!cjb node0 input
        IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Revise the physical speed of light (',CVAC,  &
                           ' in a.u.) ?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (istde,*) 'Enter the revised value:'
            READ *,C
         ELSE
            C = CVAC
         ENDIF
        ENDIF
        CALL MPI_Bcast (C, 1, MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_WORLD, ierr)
      ELSE
         C = CVAC
      ENDIF
!
! Treat some CSF's as perturbation ? Broadcasting is used only in
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      IF (NDEF .NE. 0) THEN
         IF (myid .EQ. 0) THEN
            WRITE (istde,*) 'Treat contributions of some CSFs',       &
                            ' as first-order perturbations?'
            LFORDR = GETYN ()
         ENDIF
         CALL MPI_Bcast (LFORDR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      ELSE
         LFORDR = .FALSE.
      ENDIF
 
! Get iccutblk() from the user-input
 
      IF (.NOT. LFORDR) THEN
         !...Default first
         DO i = 1, nblock
           iccutblk(i) = ncfblk(i)
         ENDDO
      ELSE
 
      ! Let master do the i/o, then broadcast
         IF (myid .EQ. 0) THEN
            WRITE (istde,*) 'There are ', nblock, 'blocks. They are:'
            WRITE (istde,*) '  block     J Parity     No of CSFs'
            DO i = 1, nblock
               WRITE (istde,*) i, idblk(i)(1:5), ncfblk(i)
            ENDDO
 
            WRITE (istde,*)
            WRITE (istde,*) 'Enter iccut for each block'
            DO jblock = 1, nblock
               WRITE (istde,*) 'Block ', jblock, '   ncf = ',        &
                              ncfblk(jblock)                         &
                      , ' id = ', idblk(jblock)(1:5)
  123          READ (istdi,*) ntmp
               IF (ntmp .GE. 0 .AND. ntmp .LE. ncfblk(jblock)) THEN
                  iccutblk(jblock) = ntmp
               ELSE
                  WRITE (istde,*) 'ICCUT out of range, re-enter:'
                  GOTO 123
               ENDIF
               write(734,*) ntmp,'! ICCUT for block',jblock
            ENDDO
         ENDIF

         CALL MPI_Bcast (iccutblk,nblock,MPI_INTEGER,0,            &
                                               MPI_COMM_WORLD,ierr)
      ENDIF ! .NOT. LFORDR
 
!*****************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN
 
!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()
 
!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.
 
!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
!*****************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Include contribution of H (Transverse)?'
         LTRANS = GETYN ()
         WRITE (istde,*) 'Modify all transverse photon frequencies?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (istde,*) 'Enter the scale factor:'
            READ *, WFACT
         ELSE
            WFACT = 1.0D00
         ENDIF
      ENDIF ! myid .EQ. 0
      CALL MPI_Bcast (LTRANS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (WFACT, 1, MPI_DOUBLE_PRECISION, 0,            &
                           MPI_COMM_WORLD, ierr)
!
! Other interactions ? One logical for each case. Done altogether
!
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Include H (Vacuum Polarisation)?'
         LVP = GETYN ()
 
         WRITE (istde,*) 'Include H (Normal Mass Shift)?'
         LNMS = GETYN ()
 
         WRITE (istde,*) 'Include H (Specific Mass Shift)?'
         LSMS = GETYN ()
 
         WRITE (istde,*) 'Estimate self-energy?'
         LSE = GETYN ()
         IF (LSE.EQV..TRUE.) THEN
             NQEDCUT = 1
            WRITE (istde,*)                                           &
        'Largest n quantum number for including self-energy for orbital'
            WRITE (istde,*) 'n should be less or equal 8'
            READ *, NQEDMAX
         ELSE
            NQEDCUT = 0
         END IF

         IF (LTRANS) THEN
          WRITE(734,'(a)')'y            ! Contribution of H Transverse?'
         ELSE
          WRITE(734,'(a)')'n            ! Contribution of H Transverse?'
         END IF
         IF (YES) THEN
           WRITE(734,'(a)') 'y            ! Modify photon frequencies?'
           WRITE(734,*) WFACT,'! Scale factor'
         ELSE
           WRITE(734,'(a)') 'n            ! Modify photon frequencies?'
         END IF
         IF (LVP) THEN
           WRITE(734,'(a)') 'y            ! Vacuum polarization?'
         ELSE
           WRITE(734,'(a)') 'n            ! Vacuum polarization?'
         END IF
         IF (LNMS) THEN
           WRITE(734,'(a)') 'y            ! Normal mass shift?'
         ELSE
           WRITE(734,'(a)') 'n            ! Normal mass shift?'
         END IF
         IF (LSMS) THEN
           WRITE(734,'(a)') 'y            ! Specific mass shift?'
         ELSE
           WRITE(734,'(a)') 'n            ! Specific mass shift?'
         END IF
         IF (LSE) THEN
           WRITE(734,'(a)') 'y            ! Self energy?'
           WRITE(734,*) NQEDMAX, '! Max n for including self energy'
         ELSE
           WRITE(734,'(a)') 'n            ! Self energy?'
         END IF
      ENDIF ! myid .EQ. 0

      CALL MPI_Bcast (LVP,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LNMS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LSMS, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (LSE,  1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
!
! Parameters controlling the radial grid
!
!cff-cjb grid
!cff use default grid
!
! 99  IF (NPARM .EQ. 0) THEN
!        RNT = EXP (-65.D0/16.D0) / Z
!        H = 0.5D0**4
!        N = MIN (220,NNNP)
!     ELSE
!CFF     .. should be Z-dependent
         RNT = 2.0D-06/Z
         H = 5.D-2
         N = NNNP
!     ENDIF
      HP = 0.D0
!     IF ( NDEF.NE.0) THEN
!        IF (myid .EQ. 0) THEN
!           WRITE (istde,*) 'The default radial grid parameters',      &
!                          ' for this case are:'
!           WRITE (istde,*) ' RNT = ',RNT,';'
!           WRITE (istde,*) ' H = ',H,';'
!           WRITE (istde,*) ' HP = ',HP,';'
!           WRITE (istde,*) ' N = ',N,';'
!           WRITE (istde,*) ' revise these values?'
!           YES = GETYN ()
!        ENDIF
! 
!         ...To prevent subsequent BCAST when YES is false
!        CALL MPI_Bcast (YES, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!        IF (YES) THEN
!           IF (myid .EQ. 0) THEN
!              WRITE (istde,*) 'Enter RNT:'
!              READ *, RNT
!              WRITE (istde,*) 'Enter H:'
!              READ *, H
!              WRITE (istde,*) 'Enter HP:'
!              READ *, HP
!              WRITE (istde,*) 'Enter N:'
!              READ *, N
!           ENDIF
!           CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (H,   1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (HP,  1, MPI_DOUBLE_PRECISION, 0,    &
!                               MPI_COMM_WORLD, ierr)
!           CALL MPI_Bcast (N,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!        ENDIF
!     ENDIF ! NDEF.NE.0
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
!   Generate $- r \times V_ (r)$
!
      CALL NUCPOT
!
!   Load the radial wavefunctions
!
      CALL SETRWFmpi (rwffile)
!
!   Write the basic parameters of the model electron cloud to the
!   .res  file; this is the second record on the file --- the
!   first is the header (see SUBROUTINE SETRES)
!
      WRITE (imcdf) NELEC, NCF, NW, nblock ! ncf is ncftot, from setcsll
!
!   Write the nuclear parameters and $- r \times V_ (r)$
!   to the  .res  file; these are the third, fourth, and fifth
!   records on the file
!
      WRITE (imcdf) Z, EMN
      WRITE (imcdf) NPARM,(PARM(I), I = 1, NPARM)
      WRITE (imcdf) N, (ZZ(I), I = 1, N), NNUC
!
!   Write the physical effects specification to the  .res  file.
!   iccutblk() is now an array of length nblock.
!   This is the sixth record on the file
!
      WRITE (imcdf) C, LFORDR, (ICCUTblk(i), i = 1, nblock),           &
                    LTRANS, WFACT, LVP, LNMS, LSMS
!
!   Write the grid data to the  .res  file; this is the seventh
!   record on the file
!
      NP10 = N + 10
      WRITE (imcdf) RNT, H, HP, (R(I), I = 1, NP10),                   &
                 (RP(I), I = 1, NP10), (RPOR(I), I = 1, NP10) ! (imcdf)
!
!   Write out the interpolated radial wavefunctions; there are
!   2*NW such records; thus the total number of records written
!   at this point is 7+2*NW
!
      DO J = 1, NW
         WRITE (imcdf) E(J), GAMA(J), PZ(J), MF(J)
         WRITE (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
      ENDDO
!
      RETURN
      END
