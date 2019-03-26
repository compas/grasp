!***********************************************************************
!                                                                      *
      SUBROUTINE SETCOF(EOL, J)
!                                                                      *
!   This  subroutine  sets  up the coefficients and orbital pointers   *
!   for the direct and exchange potentials for orbital  J .  It also   *
!   sets  up  the  coefficients  and  pointers for the inhomogeneous   *
!   terms arising from off-diagonal I (a,b) integrals.                 *
!                                                                      *
!   Call(s) to: [LIB92]: alloc, dalloc.                                *
!               [RSCF92]: alcsca, dsubrs.                              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
!   Modified by Xinghong He                 Last update: 21 Dec 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:21:02   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB
      USE memory_man
      USE orb_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE MCPA_C
      USE mpi_C
      USE pos_C
      USE scf_C
!      USE cons_C,           ONLY: EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dsubrs_I
      USE fco_I
      USE gco_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      LOGICAL  :: EOL
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      CHARACTER*6, PARAMETER :: MYNAME = 'SETCOF'
      REAL, PARAMETER :: EPS = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: INDEXS
      INTEGER :: NDIM, NAKJ, NKJJ, ILABEL, NWTERM, IB, K, NB, IR, NKJIB, KMAX, &
         KMIN, NFILE, I, JBLOCK, JBLOCKT, IENDCDUM, NCOEFF, IOS, LAB, NCONTR, &
         IA, IC, ITHIS, IRANK, IIND, IL, IORB, IYO1, IYO2, IFOUND, LOC1, LOC2, &
         ITHIS2, INDIND, NELMNTGG
      REAL(DOUBLE) :: UCFJ, SUMR, YKAB, XKAB, SUM,  CONTR
      CHARACTER :: MCPLAB*3, IDSTRING*3, MSG*128
!-----------------------------------------------
!** Locals
!     POINTER (PCOEFF,COEFF(1))
!     POINTER (PICLMN,ICLMN(1))
!     POINTER (PINDEX,indx(1))
!
!     POINTER (PNIROW,IROW(1))
!     POINTER (PNTRDA,DA(1))
!     POINTER (PNTRXA,XA(1))
!     POINTER (PNTRYA,YA(1))
!     POINTER (PNTNDA,NDA(1))
!     POINTER (PNTNXA,NXA(1))
!     POINTER (PNTNYA,NYA(1))
!
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (pnevblk, nevblk(1))
!     POINTER (pncmaxblk, ncmaxblk(1))
!
!
!     POINTER (pncfpast, ncfpast(1))
!     POINTER (pncminpast, ncminpast(1))
!     POINTER (pnevecpast, nevecpast(1))
!
      REAL(DOUBLE), DIMENSION(:), POINTER :: COEFF
      INTEGER, DIMENSION(:), pointer :: iclmn, indx
!     INTEGER, DIMENSION(:), pointer :: nda, nxa, nya
!     INTEGER, DIMENSION(:), pointer :: ncfblk, nevblk, ncmaxblk, &
!                                       ncfpast, ncminpast, nevecpast

!=======================================================================
!   Initializations
!=======================================================================

      NDIM = 1
      CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETCOF' )
      CALL ALLOC (ICLMN, NDIM, 'ICLMN', 'SETCOF')
      CALL ALLOC (INDX, NDIM, 'INDX', 'SETCOF')

      NDCOF = 0
      NXCOF = 0
      NYCOF = 0

      NAKJ = NAK(J)
      NKJJ = NKJ(J)
      UCFJ = UCF(J)

!=======================================================================
!   Generate YA coefficients that do not require MCP output list.
!   Computation distributed and then collected.
!=======================================================================

      ILABEL = 0
      NWTERM = KEY*(KEY + 1)
      DO IB = 1, NW
         ILABEL = ILABEL + NWTERM
         IF (IB == J) THEN
            KMAX = NKJJ - 1
         ELSE
            KMAX = 0
         ENDIF
         DO K = 0, KMAX, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            SUMR = 0.D0
            DO NB = 1, NBLOCK
               DO IR = MYID + 1, NCFBLK(NB), NPROCS
                  SUMR = SUMR + DSUBRS(EOL,IR,IR,NB)*FCO(K,IR + NCFPAST(NB),J,&
                     IB)
               END DO
            END DO
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            IF (IB == J) THEN
               YKAB = 2.D0*SUMR/UCFJ
            ELSE
               YKAB = SUMR/UCFJ
            ENDIF

           !*** The following IF has to be removed since YKAB
           !*** is incomplete in multiprocessor case
            !IF (ABS (YKAB) .GT. EPS) THEN
            NYCOF = NYCOF + 1
            IF (NYCOF .GT. NYDIM ) THEN
               IF (NYDIM .GT. 0) THEN
                  NYDIM = 2*NYDIM
                  CALL RALLOC(NYA, NYDIM, 'NYA', 'SETCOF')
                  CALL RALLOC(YA, NYDIM, 'YA', 'SETCOF')
               ELSE
                  NYDIM = 64
                  CALL ALLOC(NYA, NYDIM, 'NYA', 'SETCOF')
                  CALL ALLOC(YA, NYDIM, 'YA', 'SETCOF')
               ENDIF
            ENDIF
            YA(NYCOF) = YKAB
            NYA(NYCOF) = K + ILABEL
            !ENDIF

         END DO
      END DO

!=======================================================================
!   Generate XA coefficients that do not require MCP output list
!   Computation distributed and then collected.
!=======================================================================

      ILABEL = KEY*J
      NWTERM = KEY*KEY*(KEY + 1)
      DO IB = 1, NW
         ILABEL = ILABEL + NWTERM
         IF (IB == J) CYCLE
         NKJIB = NKJ(IB)
         IF (NAKJ*NAK(IB) > 0) THEN
            KMIN = ABS((NKJJ - NKJIB)/2)
         ELSE
            KMIN = ABS((NKJJ - NKJIB)/2) + 1
         ENDIF
         KMAX = (NKJJ + NKJIB)/2

         DO K = KMIN, KMAX, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            SUMR = 0.D0
            DO NB = 1, NBLOCK
               DO IR = MYID + 1, NCFBLK(NB), NPROCS
                  SUMR = SUMR + DSUBRS(EOL,IR,IR,NB)*GCO(K,IR + NCFPAST(NB),J,&
                     IB)
               END DO
            END DO
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            XKAB = SUMR/UCFJ

           !*** The following IF has to be removed since YKAB
           !*** is incomplete in multiprocessor case
            !IF (ABS (XKAB) .GT. EPS) THEN
            NXCOF = NXCOF + 1
            IF (NXCOF > NXDIM) THEN
               IF (NXDIM > 0) THEN
                  NXDIM = 2*NXDIM
                  CALL RALLOC (XA, NXDIM, 'XA', 'SETCOF')
                  CALL RALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
               ELSE
                  NXDIM = 64
                  CALL ALLOC (XA, NXDIM, 'XA', 'SETCOF')
                  CALL ALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
               ENDIF
            ENDIF
            XA(NXCOF) = XKAB
            NXA(NXCOF) = K + ILABEL
            !ENDIF
         END DO
      END DO

!=======================================================================
!   Subroutine setham called from matrix had gone through the mcp
!   files once. setcof will do it again. (Go to setmcp to see
!   contents of these records.
!=======================================================================

      DO NFILE = 30, 32 + KMAXF
         REWIND (NFILE)
         IF (NFILE /= 30) THEN
            DO I = 1, 3
               READ (NFILE)
               CYCLE
            END DO
         ELSE
            DO I = 1, 3
               READ (NFILE)
               READ (NFILE)
            END DO
         ENDIF
      END DO

!=======================================================================
!   Generate DA coefficients; these arise from the one-electron
!   integrals
!=======================================================================

      DO JBLOCK = 1, NBLOCK

         !*** Read in IROW from file mcp.30 ***
         READ (30) MCPLAB, JBLOCKT, NCF
         IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) MYNAME, '1: jblockt .NE. jblock'
            STOP
         ENDIF
         READ (30) NELMNTGG
         NELMNT = INT8(NELMNTGG)
         CALL ALLOC (IROW, NELMNT, 'IROW', 'SETCOF')
         READ (30) (IENDCDUM,I=MYID + 1,NCF,NPROCS), (IROW(I),I=1,NELMNT)

         !*** Block info file mcp.31 ***
         READ (31) MCPLAB, JBLOCKT, NCF, NCOEFF
         IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) MYNAME, '2: jblockt .NE. jblock'
            STOP
         ENDIF

         !*** Loop over labels having non-zero coefficients
         !*** it exits when no more labels for the block

         L123: DO WHILE(.TRUE.)
            READ (31, IOSTAT=IOS) LAB, NCONTR

            ! 0, 0 marks the end of a block. This is the normal exit

            IF (LAB==0 .AND. NCONTR==0) THEN
               CALL DALLOC (IROW, 'IROW', 'SETCOF')
               EXIT                              ! Actually to next block
            ENDIF

            !*** Decode the labels of I(ab) ***
            IA = MOD(LAB,KEY)
            IB = LAB/KEY

            ! At least one orbital should be J in order to have
            ! non-zero value; otherwise, goto next label.

            IF (IA/=J .AND. IB/=J) THEN
               READ (31)                         ! No contributions from this integral; skip
               CYCLE                             ! to next label
            ENDIF

            IF (NCONTR > NDIM) THEN
               CALL DALLOC (COEFF, 'COEFF', 'SETCOF')
               CALL DALLOC (ICLMN, 'ICLMN', 'SETCOF')
               CALL DALLOC (INDX, 'INDX', 'SETCOF')
               NDIM = NCONTR
               CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETCOF' )
               CALL ALLOC (ICLMN, NDIM, 'ICLMN', 'SETCOF')
               CALL ALLOC (INDX, NDIM, 'INDX', 'SETCOF')
            ENDIF

            ! Read the column index, the sparse matrix index, and the
            ! coefficient for all contributions from this integral

            READ (31) (ICLMN(I),INDX(I),COEFF(I),I=1,NCONTR)

            ! Add up all the contributions from this integral;
            ! off-diagonal contributions have double the weight

            SUM = 0.D0
            DO I = 1, NCONTR
               IR = IROW(INDX(I))
               IC = ICLMN(I)

               CONTR = DSUBRS(EOL,IR,IC,JBLOCK)*COEFF(I)

               IF (IR /= IC) CONTR = CONTR + CONTR
               SUM = SUM + CONTR
            END DO
            SUM = 0.5D0*SUM/UCFJ

            ! Put coefficients in the list. Since there is always
            ! (almost) some repeatence from different blocks, a check
            ! and merge is performed. This will significantly reduce
            ! the NDCOF and thus the number of calls to YZK later.

            !*** Find the right counting parameter ***
            IF (IA == J) THEN
               ITHIS = IB
            ELSE
               ITHIS = IA
            ENDIF

            !*** Check it against the previously recorded ***
            IF (JBLOCK > 1) THEN
               DO I = 1, NDCOF
                  IF (NDA(I) /= ITHIS) CYCLE     ! found, add the value
                  DA(I) = DA(I) + SUM
                  CYCLE  L123
               END DO
            ENDIF

            !*** Not found in the record, add an item ***
            NDCOF = NDCOF + 1
            IF (NDCOF > NDDIM) THEN
               IF (NDDIM > 0) THEN
                  NDDIM = 2*NDDIM
                  CALL RALLOC (DA, NDDIM, 'DA', 'SETCOF')
                  CALL RALLOC (NDA, NDDIM, 'NDA', 'SETCOF')
               ELSE
                  NDDIM = 64
                  CALL ALLOC (DA, NDDIM, 'DA', 'SETCOF')
                  CALL ALLOC (NDA, NDDIM, 'NDA', 'SETCOF')
               ENDIF
            ENDIF
            DA(NDCOF) = SUM                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
            NDA(NDCOF) = ITHIS

         END DO L123                             ! For labels
      END DO                                     ! For blocks

!=======================================================================
!   Generate YA and XA coefficients; these arise from the two-electron
!   integrals
!=======================================================================

      DO NFILE = 32, 32 + KMAXF

!         ...Re-position file mcp.30

         REWIND (30)
         DO I = 1, 6
            READ (30)
         END DO

!=======================================================================
!   Loop over blocks again, this time, for V-coefficients
!=======================================================================

         DO JBLOCK = 1, NBLOCK
!           ...Read in IROW from file mcp.30
            READ (30) MCPLAB, JBLOCKT, NCF
            IF (JBLOCKT /= JBLOCK) THEN
               WRITE (ISTDE, *) MYNAME, ':3 jblockt .NE. jblock'
               STOP
            ENDIF
            READ (30) NELMNTGG
            NELMNT = INT8(NELMNTGG)
            CALL ALLOC (IROW, NELMNT, 'IROW', 'SETCOF')
            READ (30) (IENDCDUM,I=MYID + 1,NCF,NPROCS), (IROW(I),I=1,NELMNT)

            READ (NFILE) MCPLAB, JBLOCKT, NCF, NCOEFF
            IF (JBLOCKT /= JBLOCK) THEN
               WRITE (ISTDE, *) MYNAME, ':4 jblockt .NE. jblock'
               STOP
            ENDIF

            K = NFILE - 32                       ! multipolarity of the integral

!=======================================================================
!   Attempt to read another block of data
!=======================================================================

  999       CONTINUE
            READ (NFILE, IOSTAT=IOS) LAB, NCONTR
!
            IF (LAB==0 .AND. NCONTR==0) THEN
               CALL DALLOC (IROW, 'IROW', 'SETCOF')
               CYCLE
            ENDIF
            !***                       k
            !*** Decode the labels of R (abcd)
            INDEXS(4) = MOD(LAB,KEY)
            LAB = LAB/KEY
            INDEXS(2) = MOD(LAB,KEY)
            LAB = LAB/KEY
            INDEXS(3) = MOD(LAB,KEY)
            INDEXS(1) = LAB/KEY

            !*** Determine the number of indices that match
            IRANK = 0
            IRANK = IRANK + COUNT(INDEXS==J)

            IF (IRANK == 0) THEN
               READ (NFILE)
               GO TO 999
            ENDIF

            !*** At least one subshell index matches; allocate storage
            !*** for reading in the rest of this block
            IF (NCONTR > NDIM) THEN
               CALL DALLOC (COEFF, 'COEFF', 'SETCOF')
               CALL DALLOC (ICLMN, 'ICLMN', 'SETCOF')
               CALL DALLOC (INDX, 'INDX', 'SETCOF')
               NDIM = NCONTR
               CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETCOF')
               CALL ALLOC (ICLMN, NDIM, 'ICLMN', 'SETCOF')
               CALL ALLOC (INDX, NDIM, 'INDX', 'SETCOF')
            ENDIF

            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            READ (NFILE) (ICLMN(I),INDX(I),COEFF(I),I=1,NCONTR)

            !*** Add up all the contributions from this integral;
            !*** off-diagonal contributions have double the weight
            SUM = 0.D0
            DO I = 1, NCONTR
               IR = IROW(INDX(I))
               IC = ICLMN(I)
               CONTR = DSUBRS(EOL,IR,IC,JBLOCK)*COEFF(I)
               IF (IR /= IC) CONTR = CONTR + CONTR
               SUM = SUM + CONTR
            END DO
            SUM = 0.5D0*SUM/UCFJ

            SELECT CASE (IRANK)
            CASE (1)

!=======================================================================
!   One matching index: exchange potential contribution
!=======================================================================

               !*** Similar to DA, find ithis ***
               ITHIS = -911                      ! initialize to an impossible value
                              ! though not necessary
               DO IIND = 1, 4
                  IF (INDEXS(IIND) /= J) CYCLE   ! at least one
                  IL = IIND + 2
                  IF (IL > 4) IL = IL - 4
                  IORB = INDEXS(IL)
                  IL = IIND + 1
                  IF (IL > 4) IL = IL - 4
                  IYO1 = INDEXS(IL)
                  IL = IIND + 3
                  IF (IL > 4) IL = IL - 4
                  IYO2 = INDEXS(IL)
                  ITHIS = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K
                  EXIT
               END DO

               IF (ITHIS == (-911)) STOP 'ithis .EQ. -911'

               !*** Check ithis against the previously recorded ***
               IF (JBLOCK > 1) THEN
                  DO I = 1, NXCOF
                     IF (NXA(I) /= ITHIS) CYCLE
                     XA(I) = XA(I) + SUM
                     GO TO 999
                  END DO
               ENDIF

               !*** Not found in records, add an item ***
               NXCOF = NXCOF + 1
               IF (NXCOF > NXDIM) THEN
                  IF (NXDIM > 0) THEN
                     NXDIM = 2*NXDIM
                     CALL RALLOC (XA, NXDIM, 'XA', 'SETCOF')
                     CALL RALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                  ELSE
                     NXDIM = 64
                     CALL ALLOC (XA, NXDIM, 'XA', 'SETCOF')
                     CALL ALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                  ENDIF
               ENDIF
               XA(NXCOF) = SUM
               NXA(NXCOF) = ITHIS

            CASE (2)

!=======================================================================
!   Two matching indices: either direct or exchange potential
!   contribution
!=======================================================================

               IFOUND = 0
               DO IIND = 1, 4
                  IF (INDEXS(IIND) /= J) CYCLE
                  IF (IFOUND == 0) THEN
                     LOC1 = IIND
                     IFOUND = IFOUND + 1
                  ELSE IF (IFOUND == 1) THEN
                     LOC2 = IIND
                     EXIT
                  ENDIF
               END DO

               IF (LOC2 - LOC1 == 2) THEN
!
!   Direct contribution
!
                  !*** Find ithis ***
                  IL = LOC1 + 3
                  IF (IL > 4) IL = IL - 4
                  IYO2 = INDEXS(IL)
                  IL = LOC1 + 1
                  IF (IL > 4) IL = IL - 4
                  IYO1 = INDEXS(IL)
                  ITHIS = (IYO2*KEY + IYO1)*KEY + K

                  !*** Check it against the previously recorded ***
                  IF (JBLOCK > 1) THEN
                     DO I = 1, NYCOF
                        IF (NYA(I) /= ITHIS) CYCLE
                        YA(I) = YA(I) + SUM + SUM
                        GO TO 999
                     END DO
                  ENDIF

                  !*** Not found, add an item ***
                  NYCOF = NYCOF + 1
                  IF (NYCOF > NYDIM) THEN
                     IF (NYDIM > 0) THEN
                        NYDIM = 2*NYDIM
                        CALL RALLOC(NYA, NYDIM, 'NYA', 'SETCOF')
                        CALL RALLOC(YA, NYDIM, 'YA', 'SETCOF')
                     ELSE
                        NYDIM = 64
                        CALL ALLOC(NYA, NYDIM, 'NYA', 'SETCOF 2b')
                        CALL ALLOC(YA, NYDIM, 'YA', 'SETCOF')
                     ENDIF
                  ENDIF
                  YA(NYCOF) = SUM + SUM
                  NYA(NYCOF) = ITHIS

               ELSE
!
!   Exchange contribution
!
                  !*** Find ithis ***
                  IL = LOC1 + 2
                  IF (IL > 4) IL = IL - 4
                  IORB = INDEXS(IL)
                  IL = LOC1 + 1
                  IF (IL > 4) IL = IL - 4
                  IYO1 = INDEXS(IL)
                  IL = LOC1 + 3
                  IF (IL > 4) IL = IL - 4
                  IYO2 = INDEXS(IL)
                  ITHIS = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K

                  !*** Check it against the previously recorded ***
                  IF (JBLOCK > 1) THEN
                     DO I = 1, NXCOF
                        IF (NXA(I) /= ITHIS) CYCLE
                        XA(I) = XA(I) + SUM + SUM
                        GO TO 999
                     END DO
                  ENDIF

                  !*** Not found, add an item ***
                  NXCOF = NXCOF + 1
                  IF (NXCOF > NXDIM) THEN
                     IF (NXDIM > 0) THEN
                        NXDIM = 2*NXDIM
                        CALL RALLOC (XA, NXDIM, 'XA', 'SETCOF')
                        CALL RALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                     ELSE
                        NXDIM = 64
                        CALL ALLOC (XA, NXDIM, 'XA', 'SETCOF')
                        CALL ALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                     ENDIF
                  ENDIF
                  XA(NXCOF) = SUM + SUM
                  NXA(NXCOF) = ITHIS

               ENDIF

            CASE (3)

!=======================================================================
!   Three matching indices: direct and exchange potential contributions
!=======================================================================

               !*** Find ithis AND ithis2
               ITHIS = -911
               ITHIS2 = -911
               DO IIND = 1, 4
                  IF (INDEXS(IIND) == J) CYCLE
                  INDIND = INDEXS(IIND)
                  IYO2 = INDIND
                  IYO1 = J
                  ITHIS = (IYO2*KEY + IYO1)*KEY + K
                  IORB = INDIND
                  IYO1 = J
                  IYO2 = J
                  ITHIS2 = ((IORB*KEY + IYO2)*KEY + IYO1)*KEY + K
               END DO

               IF (ITHIS==(-911) .OR. ITHIS2==(-911)) STOP 'ithis2'

               !*** Check the previously recorded for YA
               IF (JBLOCK > 1) THEN
                  DO I = 1, NYCOF
                     IF (NYA(I) /= ITHIS) CYCLE
                     YA(I) = YA(I) + SUM + SUM
                     GO TO 456
                  END DO
               ENDIF

               ! Not found, add an item ***
               NYCOF = NYCOF + 1
               IF (NYCOF > NYDIM) THEN
                  IF (NYDIM > 0) THEN
                     NYDIM = 2*NYDIM
                     CALL RALLOC (NYA, NYDIM, 'NYA', 'SETCOF')
                     CALL RALLOC (YA, NYDIM, 'YA', 'SETCOF')
                  ELSE
                     NYDIM = 64
                     CALL ALLOC (NYA, NYDIM, 'NYA', 'SETCOF')
                     CALL ALLOC (YA, NYDIM, 'YA', 'SETCOF')
                  ENDIF
               ENDIF
               YA(NYCOF) = SUM + SUM
               NYA(NYCOF) = ITHIS

  456          CONTINUE

               !*** Check the previously recorded for XA
               IF (JBLOCK > 1) THEN
                  DO I = 1, NXCOF
                     IF (NXA(I) /= ITHIS2) CYCLE
                     XA(I) = XA(I) + SUM
                     GO TO 999
                  END DO
               ENDIF

               ! Not found, add an item ***
               NXCOF = NXCOF + 1
               IF (NXCOF > NXDIM) THEN
                  IF (NXDIM > 0) THEN
                     NXDIM = 2*NXDIM
                     CALL RALLOC (XA, NXDIM, 'XA', 'SETCOF')
                     CALL RALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                  ELSE
                     NXDIM = 64
                     CALL ALLOC (XA, NXDIM, 'XA', 'SETCOF')
                     CALL ALLOC (NXA, NXDIM, 'NXA', 'SETCOF')
                  ENDIF
               ENDIF
               XA(NXCOF) = SUM
               NXA(NXCOF) = ITHIS2

            CASE (4)

!=======================================================================
!   Four matching indices: direct potential contribution
!=======================================================================

               !*** Find ithis AND ithis2
               IYO2 = J
               IYO1 = J
               ITHIS = (IYO2*KEY + IYO1)*KEY + K

               !*** Check the previously recorded for YA
               IF (JBLOCK > 1) THEN
                  DO I = 1, NYCOF
                     IF (NYA(I) /= ITHIS) CYCLE
                     YA(I) = YA(I) + 4.D0*SUM
                     GO TO 999
                  END DO
               ENDIF

               ! Not found, add an item ***
               NYCOF = NYCOF + 1
               IF (NYCOF > NYDIM) THEN
                  IF (NYDIM > 0) THEN
                     NYDIM = 2*NYDIM
                     CALL RALLOC(NYA, NYDIM, 'NYA', 'SETCOF')
                     CALL RALLOC(YA, NYDIM, 'YA', 'SETCOF')
                  ELSE
                     NYDIM = 64
                     CALL ALLOC(NYA, NYDIM, 'NYA', 'SETCOF 4b')
                     CALL ALLOC(YA, NYDIM, 'YA', 'SETCOF')
                  ENDIF
               ENDIF
               YA(NYCOF) = 4.D0*SUM
               NYA(NYCOF) = ITHIS

            END SELECT

            GO TO 999
         END DO                                  ! loop for V-Coefficients
      END DO

!=======================================================================
!   Deallocate storage for arrays local to this routine
!=======================================================================
      CALL DALLOC (COEFF, 'COEFF', 'SETCOF')
      CALL DALLOC (ICLMN, 'ICLMN', 'SETCOF')
      CALL DALLOC (INDX, 'INDX', 'SETCOF')

      RETURN
      END SUBROUTINE SETCOF
