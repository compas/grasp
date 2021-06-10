!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIX(dvdfirst)
!                                                                      *
!   Calls routines to form the Hamiltonian matrix and to diagonalise   *
!   it. The total angular momenta and parity of each  ASF  is found;   *
!   the eigenvectors are normalised so  that the sign of the largest   *
!   element is positive.                                               *
!                                                                      *
!   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
!               [RSCF92]: SETHAM, HMOUT.                               *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
!                                                                      *
!                                         Last revision: 24 Dec 1992   *
! Block version by Xinghong He                           07 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:52:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man
      USE damp_C
      USE def_C,           ONLY: iccmin, ncmin, ncmax
      USE DEBUG_C
      USE eigv_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE MCPA_C
      USE mpi_s
      USE orb_C
      USE pos_C
      USE peav_C
      USE syma_C,          ONLY: JPGG
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setham_I
      USE maneig_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, INTENT(IN) ::  dvdfirst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NFILE, NELEC, NTMP, JBLOCK, JBLOCKT, I, NCFPAT, NCMINPAT, &
         NEVECPAT, NCFT, NCOEFF, LAB, NCONTR, ITMP, IR, IDIAG, J, IOFSET, &
         JOTHER, IA, JOFSET, IATTMP, IASTMP, NELMNTGG
      REAL(DOUBLE), DIMENSION(:), POINTER :: CMVL
      REAL(DOUBLE) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
      INTEGER, EXTERNAL :: ddot
      LOGICAL :: FIRST
      CHARACTER :: MCPLAB*3

      SAVE FIRST
!
      DATA FIRST/ .TRUE./
!
!     POINTER (cmvl(1))
!
!-----------------------------------------------------------------------

      IF (MYID == 0) WRITE (6, *)

!   Allocate memory for CMVL once (the maximum size)
!   Save previous estimate of eigenvectors

      IF (.NOT.FIRST) THEN
         CALL ALLOC (CMVL, NVECSIZ, 'CMVL', 'MATRIX')
         CALL DCOPY (NVECSIZ, EVEC, 1, CMVL, 1)
      ENDIF

!=======================================================================
!   Position the files - MCP files (unit NFILE) for reading
!   and mixing coefficients file (unit 25) for writing
!=======================================================================

      DO NFILE = 30, 32 + KMAXF
         REWIND (NFILE)
         IF (NFILE == 30) THEN
            READ (NFILE)
            READ (NFILE)
            READ (NFILE)
         ENDIF
         READ (NFILE)
         READ (NFILE)
         READ (NFILE)
      END DO

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      IF (MYID == 0) THEN
         REWIND (25)
         READ (25)                               ! 'G92MIX'
         READ (25) NELEC, NCFTOT, NW, NTMP, NTMP, NBLOCK
         BACKSPACE (25)
         WRITE (25) NELEC, NCFTOT, NW, NCMIN, NVECSIZ, NBLOCK
      ENDIF

!=======================================================================
!   Do the job block by block
!=======================================================================

               !------------------------------------------------
      DO JBLOCK = 1, NBLOCK                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         READ (30) MCPLAB, JBLOCKT, NCF
         IF (JBLOCKT/=JBLOCK .OR. NCF/=NCFBLK(JBLOCK)) ERROR STOP &
            'matrx: jblockt .NE. jblock .OR. ncf1 .NE. ncf2'
         READ (30) NELMNTGG
         NELMNT = INT8(NELMNTGG)
         CALL ALLOC (IROW, NELMNT, 'IROW', 'MATRIX')
         CALL ALLOC (EMT, NELMNT, 'EMT', 'MATRIX')
!cjb          ALLOC (IENDC, NCF + 1,...)  ->  ALLOC (IENDC, 0, NCF,...)
!cjb     CALL ALLOC (IENDC, NCF + 1, 'IENDC', 'MATRIX' )
         CALL ALLOC (IENDC, 0, NCF, 'IENDC', 'MATRIX' )

!        ! may not be necessary if iendc is ALWAYS used
         IENDC(0:NCF) = 0                        ! the way it is assigned here.

      !...EMT will be accumulated in setham
         EMT(:NELMNT) = 0.D0

         READ (30) (IENDC(I),I=MYID + 1,NCF,NPROCS), (IROW(I),I=1,NELMNT)

         NCFPAT = NCFPAST(JBLOCK)
         NCMINPAT = NCMINPAST(JBLOCK)
         NEVECPAT = NEVECPAST(JBLOCK)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         IF (NEVBLK(JBLOCK) == 0) THEN
            DO NFILE = 31, 32 + KMAXF
               READ (NFILE) MCPLAB, JBLOCKT, NCFT, NCOEFF
               IF (JBLOCKT /= JBLOCK) ERROR STOP 'matrx: jblockt .NE. jblock'
               IF (NCFT /= NCF) ERROR STOP 'matrx: ncft .NE. ncf'

               READ (NFILE) LAB, NCONTR
               DO WHILE(LAB/=0 .OR. NCONTR/=0)
                  READ (NFILE) (ITMP,ITMP,TMP,I=1,NCONTR)
                  READ (NFILE) LAB, NCONTR
               END DO
            END DO

            CALL DALLOC (IENDC, 'IENDC', 'MATRIX')
            CALL DALLOC (IROW, 'IROW', 'MATRIX')

            CYCLE

         ENDIF

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         CALL SETHAM (JBLOCK, MYID, NPROCS)
!
!   Determine average energy
!
         EAV = 0.D0
         DO IR = myid + 1, ncf, nprocs
            EAV = EAV + EMT(IENDC(IR))
         END DO
!         DO IR = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
!            EAV = EAV + EMT(IENDC(NPROCS*(IR-1)+MYID+1))
!         END DO

         EAV = EAV/NCF
         EAVBLK(JBLOCK) = EAV

      ! Print Hamiltonian matrix and average energy
      ! hmout is not general
      !call hmout (0, 1, ncf)

         IF (MYID == 0) WRITE (*, 302) EAV

      ! Subtract the average energy from the diagonal elements
      ! to reduce the condition number of the matrix
!         DO I = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
!            EMT(IENDC(NPROCS*(I-1)+MYID+1)) = EMT(IENDC(NPROCS*(I-1)+MYID+1))&
!                - EAV
         DO i = myid + 1, ncf, nprocs
            idiag = iendc(i)     ! new mode: each row ends in diagonal
            emt(idiag) = emt(idiag) - eav
         END DO

!=======================================================================
!   Compute and store eigenpairs
!=======================================================================

         CALL MANEIG (dvdfirst, LDBPG(3),                      &
                   JBLOCK, NCFPAT, NCMINPAT, NEVECPAT, NCFTOT)

!=======================================================================
!   Damp and Schmidt orthogonalise eigenvectors for OL calculations
!=======================================================================

         IF (.NOT.FIRST) THEN

            DO J = 1, NEVBLK(JBLOCK)

               IOFSET = (J - 1)*NCF + NEVECPAT
               JOTHER = J

            ! cdamp has the original non-block feature
               CDAMPJ = CDAMP(J + NCMINPAT)
               IF (CDAMPJ == 0.D0) CYCLE         ! So SURE ???

               OMCDAJ = 1.D0 - CDAMPJ

           !...Damp eigenvector and determine the new dominant component
  123          CONTINUE
               AMAX = 0.D0
               DO I = 1, NCF
                  EVECIJ = OMCDAJ*EVEC(I+IOFSET) + CDAMPJ*CMVL(I + IOFSET)
                  EVEC(I+IOFSET) = EVECIJ
                  WA = ABS(EVECIJ)
                  IF (WA <= AMAX) CYCLE
                  AMAX = WA
                  IA = I
               END DO

            !...compute the normalization factor
               SUM = 0.D0
               DO I = 1, NCF
                  SUM = SUM + EVEC(I+IOFSET)**2
               END DO
               DNFAC = 1.D0/SQRT(SUM)

            !...Renormalize and invert as necessary
               IF (EVEC(IA+IOFSET) < 0.D0) DNFAC = -DNFAC
               CALL DSCAL (NCF, DNFAC, EVEC(IOFSET+1), 1)

            !...Schmidt orthogonalise
  234          CONTINUE
               JOTHER = JOTHER - 1
               IF (JOTHER < 1) CYCLE
               JOFSET = (JOTHER - 1)*NCF + NEVECPAT
               OVRLAP = DDOT(NCF,EVEC(IOFSET+1),1,EVEC(JOFSET+1),1)
               IF (OVRLAP /= 0.D0) THEN          ! So SURE ???
                  OMCDAJ = 1.D0
                  CDAMPJ = -OVRLAP
                  CALL DCOPY (NCF, EVEC(JOFSET+1), 1, CMVL(IOFSET + 1), 1)
                  GO TO 123
               ELSE
                  GO TO 234
               ENDIF
            END DO
         ENDIF

!   Write out the eigenpair information: ASF symmetries, eigenvalues,
!   and eigenvectors to GRASP92 mixing coefficients File

         IF (NEVBLK(JBLOCK) == 0) THEN
            IATTMP = 999
            IASTMP = 999
         ELSE
!GGGG
            iattmp = IABS(JPGG(jblock))
            IF(JPGG(jblock) >= 0) THEN
!cjb iastmp error corrected by Pawel Syty
!PS            iastmp = 2
               iastmp = 1
            ELSE
!cjb iastmp error corrected by Pawel Syty
!PS            iastmp = 1
               iastmp = -1
            END IF
!GG            IATTMP = IATJPO(NCMINPAT + 1)
!GG            IASTMP = IASPAR(NCMINPAT + 1)
         ENDIF

         IF (MYID == 0) THEN
            WRITE (25) JBLOCK, NCF, NEVBLK(JBLOCK), IATTMP, IASTMP
            WRITE (25) (ICCMIN(I + NCMINPAT),I=1,NEVBLK(JBLOCK))
            WRITE (25) EAV, (EVAL(I + NCMINPAT),I=1,NEVBLK(JBLOCK))
            WRITE (25) ((EVEC(I+(J-1)*NCF+NEVECPAT),I=1,NCF),J=1,NEVBLK(JBLOCK)&
               )
         ENDIF

         CALL DALLOC (EMT, 'EMT', 'MATRIX')
         CALL DALLOC (IENDC, 'IENDC', 'MATRIX')
         CALL DALLOC (IROW, 'IROW', 'MATRIX')

!----------------------
      END DO
!----------------------
!
!   Deallocate the temporary storage
!
      IF (.NOT.FIRST) THEN
         CALL DALLOC (CMVL, 'CMVL', 'MATRIX')
      ELSE
         FIRST = .FALSE.
      ENDIF

  302 FORMAT(' Average energy = ',1P,D18.10,' Hartrees')

      RETURN
      END SUBROUTINE MATRIX
