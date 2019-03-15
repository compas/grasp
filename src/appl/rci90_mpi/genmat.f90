!***********************************************************************
!                                                                      *
      SUBROUTINE GENMAT(ATWINV, JBLOCK, MYID, NPROCS, ELSTO, IRESTART, SLF_EN)
!
!        Generate Hamiltonian matrix for all blocks
!        This routine calls setham to do the computation. It makes
!        sure that the hamiltonian matrix is complete.
!  Only node-0 has the correct, non-zero elsto. And in any case, elsto
!  is not added to the matrix elements in files *.res
!  This routine has been re-written to work in both single- and
!  multi- processors.
!
! Xinghong He 1998-06-23
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE eigv_C
      USE iccu_C
      USE orb_C
      USE where_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE posfile_I
      USE setham_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JBLOCK
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      INTEGER, INTENT(OUT) :: IRESTART
      REAL(DOUBLE)  :: ATWINV
      REAL(DOUBLE)  :: ELSTO
      REAL(DOUBLE), DIMENSION(*) :: SLF_EN
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IREAD, IOS, NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM, I, &
         IOS2, NELC, IR, IROWDUM, IPOS, NROWS, J, ICSTRT
      REAL(DOUBLE) :: STOEL, DUM, EAV0
!-----------------------------------------------
!
      NELMNT = 0                      ! Counting continues in setham
      EAV    = 0.D0
      ELSTO  = 0.D0

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      IREAD = 0                       ! # of rows read, initialization necessary

      READ (IMCDF, IOSTAT=IOS) NCFDUM, ICCUTDUM, MYIDDUM, NPROCSDUM
      IRESTART = IOS

      IF (IOS == 0) THEN

         IF (NCF/=NCFDUM .OR. ICCUT(1)/=ICCUTDUM .OR. MYID/=MYIDDUM .OR. &
            NPROCS/= NPROCSDUM) THEN
            WRITE (ISTDE, *) NCF, NCFDUM, ICCUT(1), ICCUTDUM, MYID, MYIDDUM, &
               NPROCS, NPROCSDUM, 'check'
            STOP 'genmat:1'
         ENDIF

         DO I = MYID + 1, NCF, NPROCS
            READ (IMCDF, IOSTAT=IOS2) NELC, STOEL, (DUM,IR=2,NELC), EAV0, (&
               IROWDUM,IR=1,NELC)
                            ! Lower triangle row-mode, diagonal last
            IF (IOS2 == 0) THEN
               IREAD = IREAD + 1
               NELMNT = NELMNT + NELC
               EAV = EAV + EAV0
               ELSTO = STOEL
            ELSE
               EXIT
            ENDIF
         END DO
         IPOS = 7 + NW + NW + IREAD + 1
      ELSE
         IPOS = 7 + NW + NW
      ENDIF

! Find the maximum number of rows

      NROWS = (NCF - MYID - 1 + NPROCS)/NPROCS
      IF (NCF < NPROCS) NROWS = NCF/(MYID + 1)

! Report the number of rows read.
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

      WRITE (ISTDE, *) IREAD, ' (total ', NROWS, ') rows read from .res'
      IF (MYID == 0) WRITE (24, *) IREAD, ' (total ', NROWS, &
         ') rows read from .res'

! Position the file for the next record from setham

      DO I = 1, JBLOCK - 1
         J = (NCFBLK(I) - MYID - 1 + NPROCS)/NPROCS
         IF (NCFBLK(I) < NPROCS) J = NCFBLK(I)/(MYID + 1)
         IPOS = IPOS + J + 1
      END DO
      CALL POSFILE (0, IMCDF, IPOS)

      IF (IOS /= 0) WRITE (IMCDF) NCF, ICCUT(1), MYID, NPROCS

      IF (IREAD < NROWS) THEN
         ICSTRT = IREAD*NPROCS + MYID + 1
!     ...Generate the rest of the Hamiltonian matrix
         CALL SETHAM (MYID, NPROCS, JBLOCK, ELSTO, ICSTRT, NELMNT, ATWINV, &
            SLF_EN)
      ELSE
         NELMNTTMP = NELMNT
         NCFTMP = NCF
      ENDIF

      RETURN
      END SUBROUTINE GENMAT
