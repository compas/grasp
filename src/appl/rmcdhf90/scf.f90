!***********************************************************************
!                                                                      *
      SUBROUTINE SCF(EOL, RWFFILE2) 
!                                                                      *
!   This  subroutine  performs  the SCF iterations. The procedure is   *
!   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
!   (1986) 290.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!               [RSCF92]: improv, matrix, MAXARR, newco,               *
!                         ORBOUT, ORTHSC, setlag.                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   Block version by Xinghong He          Last revision: 05 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:16:00   1/ 5/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE memory_man
      USE blkidx_C
      USE default_C
      USE def_C
      USE debug_C
      USE eigv_C
      USE fixd_C
      USE hblock_C
      USE iounit_C
      USE lagr_C
      USE MCPA_C 
      USE mpi_s
      USE pos_c
      USE peav_C
      USE ORB_C 
      USE orba_C
      USE SCF_C
      USE ORTHCT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE matrix_I 
      USE newco_I 
      USE setlag_I 
      USE improv_I 
      USE maxarr_I 
      USE prwf_I 
      USE orthsc_I 
      USE orbout_I 
      USE endsum_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL  :: EOL 
      CHARACTER  :: RWFFILE2*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, NIT, JSEQ, KOUNT, K 
      REAL(DOUBLE) :: WTAEV, WTAEV0, DAMPMX 
      LOGICAL :: CONVG, LSORT, dvdfirst
!-----------------------------------------------


      NCFTOT = NCF 
      !IF (myid .EQ. 0) PRINT *, '===SCF==='
 
!=======================================================================
!   Determine Orthonomalization order --- lsort
!=======================================================================
 
      IF (NDEF == 0) THEN 
         LSORT = .FALSE. 
      ELSE 
  123    CONTINUE 
         WRITE (ISTDE, *) 'Orthonomalization order? ' 
         WRITE (ISTDE, *) '     1 -- Update order' 
         WRITE (ISTDE, *) '     2 -- Self consistency connected' 
         READ (ISTDI, *) J 
         IF (J == 1) THEN 
            LSORT = .FALSE. 
         ELSE IF (J == 2) THEN 
            LSORT = .TRUE. 
         ELSE 
            WRITE (ISTDE, *) 'Input is wrong, redo...' 
            GO TO 123 
         ENDIF 
      ENDIF 
 
!=======================================================================
!   Deallocate storage that will no longer be used
!=======================================================================
 
!GG      CALL DALLOC (JQSA, 'JQSA', 'SCF') 
 
!=======================================================================
!   Allocate and fill in auxiliary arrays
!=======================================================================
 
      CALL ALLOC (NCFPAST, NBLOCK, 'NCFPAST', 'SCF') 
      CALL ALLOC (NCMINPAST, NBLOCK, 'NCMINPAST', 'SCF') 
      CALL ALLOC (NEVECPAST, NBLOCK, 'NEVECPAST', 'SCF') 
      CALL ALLOC (EAVBLK, NBLOCK, 'EAVBLK', 'SCF') 
 
      NCFPAST(1) = 0 
      NCMINPAST(1) = 0 
      NEVECPAST(1) = 0 
      DO I = 2, NBLOCK 
         NCFPAST(I) = NCFPAST(I-1) + NCFBLK(I - 1) 
         NCMINPAST(I) = NCMINPAST(I-1) + NEVBLK(I - 1) 
         NEVECPAST(I) = NEVECPAST(I-1) + NEVBLK(I - 1)*NCFBLK(I - 1) 
      END DO 
 
      !*** Size of the eigenvector array for all blocks
      NVECSIZ = NEVECPAST(NBLOCK) + NEVBLK(NBLOCK)*NCFBLK(NBLOCK) 
 
      IF (EOL) THEN 
         CALL ALLOC (EVAL, NCMIN, 'EVAL', 'SCF') 
         CALL ALLOC (EVEC, NVECSIZ, 'EVEC', 'SCF') 
!GG         CALL ALLOC (IATJPO, NCMIN, 'IATJPO', 'SCF') 
!GG         CALL ALLOC (IASPAR, NCMIN, 'IASPAR', 'SCF') 
      ENDIF 
 
!=======================================================================
!
!=======================================================================
      NDDIM = 0 
      NXDIM = 0 
      NYDIM = 0 
 
!     This call should only be made AFTER the call to newco
!     CALL setlag (EOL)
 
!   For (E)OL calculations, determine the level energies and
!   mixing coefficients

!CFF   .. set the logical variable dvdfirst
      dvdfirst = .true. 
      IF (EOL) THEN 
         CALL MATRIX (dvdfirst)
         CALL NEWCO (WTAEV) 
      ENDIF 
      WTAEV0 = 0.0 
      dvdfirst = .false.
      DO NIT = 1, NSCF 
         IF (MYID == 0) WRITE (*, 301) NIT 
 
!   For all pairs constrained through a Lagrange multiplier, compute
!   the Lagrange multiplier
 
         CALL SETLAG (EOL) 
 
!   Improve all orbitals in turn
 
         DAMPMX = 0.0 
         IF (MYID == 0) WRITE (*, 302) 
         DO J = 1, NW 
            JSEQ = IORDER(J) 
            IF (LFIX(JSEQ)) CYCLE  
            CALL IMPROV (EOL, JSEQ, LSORT, DAMPMX) 
         END DO 
!
!   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
!   improve it
!
!         write(istde,*) 'nsic=',nsic
         DO KOUNT = 1, NSIC 
            CALL MAXARR (K) 
            IF (K == 0) THEN 
               CONVG = .TRUE. 
               GO TO 3 
            ELSE 
               IF (SCNSTY(K) <= ACCY) THEN 
                  CONVG = .TRUE. 
                  GO TO 3 
               ENDIF 
            ENDIF 
            CALL IMPROV (EOL, K, LSORT, DAMPMX) 
         END DO 
 
         CALL MAXARR (K) 
 
         IF (K == 0) THEN 
            CONVG = .TRUE. 
         ELSE 
            IF (SCNSTY(K) <= ACCY) THEN 
               CONVG = .TRUE. 
            ELSE 
               CONVG = .FALSE. 
            ENDIF 
         ENDIF 
 
    3    CONTINUE 
         IF (LDBPR(24) .AND. MYID==0) CALL PRWF (0) 
 
!   Perform Gram-Schmidt process
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         IF (.NOT.ORTHST) CALL ORTHSC 
 
!   Write the subshell radial wavefunctions to the .rwf file
 
         IF (MYID == 0) CALL ORBOUT (RWFFILE2) 
 
         IF (EOL) THEN 
            CALL MATRIX(dvdfirst)
            CALL NEWCO (WTAEV) 
         ENDIF 
!        Make this a relative convergence test
!        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
!    &   DAMPMX.LT.1.0D-4) CONVG=.true.
!        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
!    &             ABS((WTAEV-WTAEV0)/WTAEV)
!cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
!cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
!cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE. 
         WTAEV0 = WTAEV 
         IF (.NOT.CONVG) CYCLE  
         IF (LDBPR(25) .AND. .NOT.LDBPR(24) .AND. MYID==0) CALL PRWF (0) 
            !IF (EOL) CALL matrix (dvdfirst)
         GO TO 5 
 
      END DO 
 
      IF (MYID == 0) WRITE (ISTDE, *) ' Maximum iterations in SCF Exceeded.' 
 
    5 CONTINUE 
      DO I = 31, 32 + KMAXF 
         CLOSE(I)                                ! The MCP coefficient files 
      END DO 
 
      IF (MYID == 0) THEN 
         !CLOSE (23)     ! The .rwf file
         CLOSE(25)                               ! The .mix file 
!
!   Complete the summary - moved from rscf92 for easier alloc/dalloc
!
         CALL ENDSUM 
      ENDIF 
!
!   Deallocate storage
!
      CALL DALLOC (WT, 'WT', 'SCF')                       !Either getold or getald 
 
      IF (NEC > 0) THEN 
         CALL DALLOC (IECC, 'IECC', 'SCF') 
         CALL DALLOC (ECV, 'ECV', 'SCF') 
!GG         CALL DALLOC (IQAR, 'IQAR', 'SCF') 
         CALL DALLOC (IQA, 'IQA', 'SCF') 
      ENDIF 
 
      IF (NDDIM > 0) THEN 
         CALL DALLOC (DA, 'DA', 'SCF') 
         CALL DALLOC (NDA, 'NDA', 'SCF') 
      ENDIF 
 
      IF (NXDIM > 0) THEN 
         CALL DALLOC (XA, 'XA', 'SCF') 
         CALL DALLOC (NXA, 'NXA', 'SCF') 
      ENDIF 
 
      IF (NYDIM > 0) THEN 
         CALL DALLOC (YA, 'YA', 'SCF') 
         CALL DALLOC (NYA, 'NYA', 'SCF') 
      ENDIF 
 
      IF (EOL) THEN 
         CALL DALLOC (EVAL, 'EVAL', 'SCF') 
         CALL DALLOC (EVEC, 'EvEC', 'SCF') 
!GG         CALL DALLOC (IATJPO, 'IATJPO', 'SCF') 
!GG         CALL DALLOC (IASPAR, 'IASPAR', 'SCF') 
         CALL DALLOC (NCMAXBLK, 'NCMAXBLK', 'SCF') ! getold.f 
         CALL DALLOC (EAVBLK, 'EAVBLK', 'SCF')     ! getold.f 
         CALL DALLOC (IDXBLK, 'IDXBLK', 'SCF')     ! Allocated in getold.f 
         CALL DALLOC (ICCMIN, 'ICCMIN', 'SCF')       ! Allocated in items.f<-getold.f 
      ENDIF 
!
      CALL DALLOC (NCFPAST, 'NCFPAST', 'SCF') 
      CALL DALLOC (NCMINPAST, 'NCMINPAST', 'SCF') 
      CALL DALLOC (NEVECPAST, 'NEVECPAST', 'SCF') 
 
  301 FORMAT(/,' Iteration number ',1I3,/,' --------------------') 
  302 FORMAT(41X,'Self-            Damping'/,&
         'Subshell    Energy    Method   P0    ',&
         'consistency  Norm-1  factor  JP',' MTP INV NNP'/) 
 
      RETURN  
      END SUBROUTINE SCF 
