!***********************************************************************
!                                                                      *
      SUBROUTINE GETOLDmpi(IDBLK)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:25:01   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE memory_man
      USE blkidx_C
      USE damp_C, ONLY: cdamp, odamp
      USE def_C, rwtdum=>wt
      USE default_C
      USE fixd_C, ONLY: lfix, nfix
      USE hblock_C
      USE iounit_C
      USE invt_C, ONLY: noinvt
      USE orthct_C
      USE ORB_C
      USE ORBA_C, ONLY: IORDER
      USE CORRE_C, ONLY: LCORRE
      USE scf_C, ONLY: SCNSTY,METHOD
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lodstate_I
      USE getoldwt_I
      USE prtrsl_I
      USE getrsl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM
      INTEGER , DIMENSION(NNNW) :: INDX
      INTEGER :: I, NSUBS, NORDER, LOC, NOFFSET, JBLOCK, J
      LOGICAL :: GETYN, YES
      CHARACTER :: RECORD*256
!-----------------------------------------------

      !WRITE (istde,*) 'EOL type calculation;'

! lodstate fills
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)

      CALL ALLOC (NCMAXBLK, NBLOCK, 'NCMAXBLK', 'GETOLDmpi')
      CALL ALLOC (NEVBLK, NBLOCK, 'NEVBLK', 'GETOLDmpi' )

!cjb  ncmaxblk & nevblk initialised in lodstate
!     NEVBLK = 0
!     NCMAXBLK = 0
!cjb

!cjb LODSTATE(NBLOCK, NCFBLK, IDBLK, NEVBLK, NCMAXBLK) -> (IDBLK)
!cjb   CALL LODSTATE (NBLOCK, NCFBLK(1), IDBLK, NEVBLK, NCMAXBLK)
      IF (myid == 0) THEN
         CALL LODSTATE (IDBLK)
      END IF
!cjb

      CALL MPI_Bcast (ncmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (nevblk(1), nblock, MPI_INTEGER, 0,            &
                      MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ncmaxblk(1), nblock, MPI_INTEGER, 0,          &
                      MPI_COMM_WORLD, ierr)

      ! actually, ncmin will always be positive
      IF (ncmin /= 0) THEN
         IF (myid /=  0) THEN   ! lodstate allocated it from node-0
!GG            CALL alloc (pccmin, ncmin, 'PCCMIN', 'GETOLD')
            CALL alloc (iccmin, ncmin, 'ICCMIN', 'GETOLDmpi')
         ENDIF

         CALL MPI_Bcast (iccmin(1), ncmin, MPI_INTEGER, 0,         &
                         MPI_COMM_WORLD, ierr)
      ENDIF

!
!   Allocate the storage for and set the weights
!
      CALL ALLOC (WEIGHT, NCMIN, 'WEIGHT', 'GETOLDmpi')

      IF (myid == 0) THEN
         CALL GETOLDWT (NDEF, NCMIN, WEIGHT)
      END IF
      CALL MPI_Bcast (weight(1), ncmin, MPI_DOUBLE_PRECISION, 0,  &
                      MPI_COMM_WORLD, ierr)
!
!   Eigenvector damping
!
      CALL ALLOC (CDAMP, NCMIN, 'CDAMP', 'GETOLDmpi')
!
      CDAMP(:NCMIN) = 0.D0
!
!   Print the list of all subshells
!
      IF (myid == 0) THEN
         WRITE (ISTDE, *) 'Radial functions'
         CALL PRTRSL
      END IF
!
!   Determine which orbitals are to be varied, which are fixed.
!   Quantities determined here:
!     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
!   Instead of broadcasting these quantities, we broadcast
!   the intermediate result from GETRSL (see below)
!
      LFIX(:NW) = .TRUE.

      IF (myid == 0) THEN
         WRITE (ISTDE, *) 'Enter orbitals to be varied (Updating order)'
         CALL GETRSL (INDX, NSUBS)
      END IF
      CALL MPI_Bcast (nsubs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      IF (nsubs > 0) CALL MPI_Bcast (indx(1), nsubs, MPI_INTEGER, 0,  &
                                     MPI_COMM_WORLD, ierr)

      LFIX(INDX(:NSUBS)) = .FALSE.
!XHH      give a big value, rather than zero to scnsty()
      SCNSTY(INDX(:NSUBS)) = 1.D20
      NFIX = NW - NSUBS
      IF (NFIX == NW) THEN
         IF(MYID == 0)                                                 &
         WRITE (ISTDE,*)'All subshell radial wavefunctions are fixed;',&
                        ' performing CI calculations with RCI.'
      END IF

!   Determine orbital updating order

      NORDER = 0
      DO I = 1, NW
         IORDER(I) = I
         IF (LFIX(I)) CYCLE
         NORDER = NORDER + 1
         IORDER(I) = INDX(NORDER)
      END DO
!
!XHH added a array to store the index of the correlation functions
!
      LCORRE(:NW) = .TRUE.

      IF (myid == 0) THEN
         WRITE (ISTDE, *) 'Which of these are spectroscopic orbitals?'
         CALL GETRSL (INDX, NSUBS)
      END IF
      CALL MPI_Bcast (nsubs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      IF (NSUBS .GT. 0) THEN
         CALL MPI_Bcast (indx(1), nsubs, MPI_INTEGER, 0,            &
                         MPI_COMM_WORLD, ierr)
!GG      IF (NSUBS > 0) THEN
         DO I = 1, NSUBS
            LOC = INDX(I)
            IF (LFIX(LOC)) CYCLE
            METHOD(LOC) = 1
            NOINVT(LOC) = .FALSE.
            ODAMP(LOC) = 0.D0
            LCORRE(LOC) = .FALSE.
         END DO
      ENDIF

! Set NSIC. It will be non-zero if all orbitals to be varied are
! spectroscopic orbitals

      NSIC = 4 + (NW - NFIX)/4
      DO I = 1, NW
         IF (.NOT.(.NOT.LFIX(I) .AND. LCORRE(I))) CYCLE
         NSIC = 0
         EXIT
      END DO
!
      NSCF = 24
      NSOLV = 3
      ORTHST = .TRUE.
!
!   Make the allocation for the auxiliary vector required
!   by SUBROUTINE NEWCO
!
      CALL ALLOC (RWTDUM, NCMIN, 'RWTDUM', 'GETOLDmpi')
!
!   Place the block numbers of the all ncmin eigenstate(wanted)
!   in array idxblk
!
      CALL ALLOC (IDXBLK, NCMIN, 'IDXBLK', 'GETOLDmpi')
      NOFFSET = 0
      DO JBLOCK = 1, NBLOCK
         DO J = 1, NEVBLK(JBLOCK)
            IDXBLK(J+NOFFSET) = JBLOCK
         END DO
         NOFFSET = NOFFSET + NEVBLK(JBLOCK)
      END DO
      IF (NOFFSET /= NCMIN) ERROR STOP 'getold: ncmin trouble'

      RETURN
      END SUBROUTINE GETOLDmpi
