!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIXmpi(dvdfirst)
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
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:52:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE, LONG
      USE memory_man
      USE damp_C
      USE def_C, ONLY: iccmin, ncmin, ncmax
      USE DEBUG_C
      USE eigv_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE MCPA_C
      USE mpi_C
      USE orb_C
      USE pos_C
      USE peav_C
      USE syma_C,          ONLY: JPGG
      USE rmcdhf_mem_C,    ONLY: IMIN1_30, IMAX1_30,IMIN2_30, IMAX2_30
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setham_I
      USE maneigmpi_I
      USE read1_mem_I
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
      INTEGER(LONG), DIMENSION(20) :: NCONTR_tot, NCONTR_tot2

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
         CALL ALLOC (CMVL, NVECSIZ, 'CMVL', 'MATRIXmpi')
         CALL DCOPY (NVECSIZ, EVEC, 1, CMVL, 1)
      ENDIF

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

      NCONTR_tot = 0
      NCONTR_tot2 = 0
               !------------------------------------------------
      DO JBLOCK = 1, NBLOCK                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         NELMNTGG = IMAX2_30(JBLOCK) - IMIN2_30(JBLOCK)
         NELMNT = INT8(NELMNTGG)
         NCF = IMAX1_30(JBLOCK) - IMIN1_30(JBLOCK)
         CALL ALLOC (IROW, NELMNT, 'IROW', 'MATRIXmpi')
         CALL ALLOC (EMT, NELMNT, 'EMT', 'MATRIXmpi')
         CALL ALLOC (IENDC, NCF + 1, 'IENDC', 'MATRIXMPI' )

!                                                ! may not be necessary if iendc is ALWAYS used
!GG         IENDC(:NCF) = 0                        ! the way it is assigned here.
         IENDC(0:NCF) = 0

      !...EMT will be accumulated in setham
         EMT(:NELMNT) = 0.D0
         CALL read03_mem(JBLOCK)

         NCFPAT = NCFPAST(JBLOCK)
         NCMINPAT = NCMINPAST(JBLOCK)
         NEVECPAT = NEVECPAST(JBLOCK)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         IF (NEVBLK(JBLOCK) == 0) THEN
            DO NFILE = 31, 32 + KMAXF
               NCONTR_tot2(NFILE-30) = NCONTR_tot2(NFILE-30) + 1
               CALL read1_mem(NFILE, NCONTR_tot2(NFILE-30), LAB, NCONTR)
               DO WHILE(LAB/=0 .OR. NCONTR/=0)
                  NCONTR_tot(NFILE-30) = NCONTR_tot(NFILE-30) + NCONTR
                  NCONTR_tot2(NFILE-30) = NCONTR_tot2(NFILE-30) + 1
                  CALL read1_mem(NFILE, NCONTR_tot2(NFILE-30), LAB, NCONTR)
               END DO
            END DO

            CALL DALLOC (IENDC, 'IENDC', 'MATRIXmpi')
            CALL DALLOC (IROW, 'IROW', 'MATRIXmpi')

            CYCLE

         ENDIF

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         CALL SETHAM (JBLOCK, MYID, NPROCS, NCONTR_tot, NCONTR_tot2)
!
!   Determine average energy
!
         EAV = 0.D0
         DO IR = myid + 1, ncf, nprocs
            EAV = EAV + EMT(IENDC(IR))
         END DO
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_Allreduce (eav,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                             MPI_COMM_WORLD,ierr)
         eav = tmp / ncf
!GG         EAV = EAV/NCF
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
         CALL MANEIGmpi (dvdfirst, LDBPG(3),                      &
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
         IF(JPGG(jblock) .GE. 0) THEN
            iastmp = 1
         ELSE
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

         CALL DALLOC (EMT, 'EMT', 'MATRIXmpi')
         CALL DALLOC (IENDC, 'IENDC', 'MATRIXmpi')
         CALL DALLOC (IROW, 'IROW', 'MATRIXmpi')

!----------------------
      END DO
!----------------------
!
!   Deallocate the temporary storage
!
      IF (.NOT.FIRST) THEN
         CALL DALLOC (CMVL, 'CMVL', 'MATRIXmpi')
      ELSE
         FIRST = .FALSE.
      ENDIF

  302 FORMAT(' Average energy = ',1P,D18.10,' Hartrees')

      RETURN
      END SUBROUTINE MATRIXmpi
