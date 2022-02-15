!***********************************************************************
!                                                                      *
      SUBROUTINE GETOLD(IDBLK)
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:25:01   1/ 6/07
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

      CALL ALLOC (NCMAXBLK, NBLOCK, 'NCMAXBLK', 'GETOLD')
      CALL ALLOC (NEVBLK, NBLOCK, 'NEVBLK', 'GETOLD' )

!cjb  ncmaxblk & nevblk initialised in lodstate
!     NEVBLK = 0
!     NCMAXBLK = 0
!cjb

!cjb LODSTATE(NBLOCK, NCFBLK, IDBLK, NEVBLK, NCMAXBLK) -> (IDBLK)
!     CALL LODSTATE (NBLOCK, NCFBLK(1), IDBLK, NEVBLK, NCMAXBLK)
      CALL LODSTATE (IDBLK)
!cjb
!
!   Allocate the storage for and set the weights
!
      CALL ALLOC (WEIGHT, NCMIN, 'WEIGHT', 'GETOLD')

      CALL GETOLDWT (NDEF, NCMIN, WEIGHT)
!
!   Eigenvector damping
!
      CALL ALLOC (CDAMP, NCMIN, 'CDAMP', 'GETOLD')
!
      CDAMP(:NCMIN) = 0.D0
!
!   Print the list of all subshells
!
      WRITE (ISTDE, *) 'Radial functions'
      CALL PRTRSL
!
!   Determine which orbitals are to be varied, which are fixed.
!   Quantities determined here:
!     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
!   Instead of broadcasting these quantities, we broadcast
!   the intermediate result from GETRSL (see below)
!
      LFIX(:NW) = .TRUE.

      WRITE (ISTDE, *) 'Enter orbitals to be varied (Updating order)'
      CALL GETRSL (INDX, NSUBS)

      LFIX(INDX(:NSUBS)) = .FALSE.
!XHH      give a big value, rather than zero to scnsty()
      SCNSTY(INDX(:NSUBS)) = 1.D20
      NFIX = NW - NSUBS
      IF (NFIX == NW) WRITE (ISTDE, *) &
         'All subshell radial wavefunctions are fixed;', &
         ' performing CI calculation.'

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

      WRITE (ISTDE, *) 'Which of these are spectroscopic orbitals?'
      CALL GETRSL (INDX, NSUBS)
      IF (NSUBS > 0) THEN
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

      NSIC = (NW - NFIX)/4
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
      CALL ALLOC (RWTDUM, NCMIN, 'RWTDUM', 'GETOLD')
!
!   Place the block numbers of the all ncmin eigenstate(wanted)
!   in array idxblk
!
      CALL ALLOC (IDXBLK, NCMIN, 'IDXBLK', 'GETOLD')
      NOFFSET = 0
      DO JBLOCK = 1, NBLOCK
         DO J = 1, NEVBLK(JBLOCK)
            IDXBLK(J+NOFFSET) = JBLOCK
         END DO
         NOFFSET = NOFFSET + NEVBLK(JBLOCK)
      END DO
      IF (NOFFSET /= NCMIN) ERROR STOP 'getold: ncmin trouble'

      RETURN
      END SUBROUTINE GETOLD
