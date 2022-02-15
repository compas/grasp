!***********************************************************************
!                                                                      *
      SUBROUTINE SETHAM(JBLOCK, MYID, NPROCS, NCONTR_tot, NCONTR_tot2)
!                                                                      *
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:23:52   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE, LONG
      USE parameter_def,    ONLY: KEYORB
      USE memory_man
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE hmat_C
      USE MCPA_C
      USE orb_C
      USE pos_C
      USE iounit_C
      USE rmcdhf_mem_C,     ONLY: INDX, COEFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      USE rinti_I
      USE fco_I
      USE slater_I
      USE gco_I
      USE read1_mem_I
      USE read2_mem_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JBLOCK
      INTEGER, INTENT(IN) :: MYID, NPROCS
      INTEGER(LONG), DIMENSION(20) :: NCONTR_tot, NCONTR_tot2
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      CHARACTER*6, PARAMETER :: MYNAME = 'SETHAM'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCFPAT, IA, IR, ITMP, IDIAG, KM, K, K0, IB, NKJIA, NKJIB, KMIN&
         , KMAX, NDIM, JBLOCKT, NCFT, NCOEFF, IOS, LAB, NCONTR, ICLMNDUM, I, &
         LOC, NFILE, ID, IC
      INTEGER(LONG) :: LOD      !  INTEGER*8 version of LOC
      REAL(DOUBLE) :: DIAA, COEF, F0AA, FKAA, F0AB, GKAB, TEGRAL
      LOGICAL :: SET
      CHARACTER :: MCPLAB*3
!-----------------------------------------------
!
!     !*** needs EMT, IENDC.
!     POINTER (PNTEMT,EMT(1))
!     POINTER (PIENDC,IENDC(0:*))
!     POINTER (PNIROW,IROW(1))
!
!     !*** needs ncf, nw
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     POINTER (pncfpast, ncfpast(1))
!     POINTER (pncminpast, ncminpast(1))
!     POINTER (pnevecpast, nevecpast(1))
!
!     POINTER (PCOEFF,COEFF(1))
!     !POINTER (PICLMN,ICLMNdum)
!     POINTER (PINDX,INDX(1))
!
!-----------------------------------------------------------------------
      NCFPAT = NCFPAST(JBLOCK)
!=======================================================================
!   Accumulate diagonal terms that do not require MCP coefficients
!=======================================================================

!=======================================================================
!   Piece involving I(a,a) integrals
!=======================================================================

      DO IA = 1, NW
         SET = .FALSE.
         DO IR = MYID + 1, NCF, NPROCS
            ITMP = IQ(IA,IR + NCFPAT)
            IF (ITMP <= 0) CYCLE
            !*** Occupation number not zero ...
            IF (.NOT.SET) THEN
               DIAA = RINTI(IA,IA,0)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)                    ! lower-triangle-by-rows mode
            EMT(IDIAG) = EMT(IDIAG) + ITMP*DIAA
         END DO
      END DO

!=======================================================================
!                    0
!   Piece involving F (a,a) integrals
!=======================================================================

      DO IA = 1, NW
         SET = .FALSE.
         DO IR = MYID + 1, NCF, NPROCS
            COEF = FCO(0,IR + NCFPAT,IA,IA)
            IF (COEF == 0.D0) CYCLE
            !*** Angular coefficient not zero ...
            IF (.NOT.SET) THEN
               F0AA = SLATER(IA,IA,IA,IA,0)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)
            EMT(IDIAG) = EMT(IDIAG) + COEF*F0AA
         END DO
      END DO

!=======================================================================
!                    k
!   Piece involving F (a,a) integrals
!=======================================================================

      KM = 0
      K = 0
    6 CONTINUE
      K = K + 2
      DO IA = 1, NW
         K0 = NKJ(IA) - 1
         KM = MAX0(K0,KM)
         IF (K > K0) CYCLE
         SET = .FALSE.

         DO IR = MYID + 1, NCF, NPROCS
            COEF = FCO(K,IR + NCFPAT,IA,IA)
            IF (COEF == 0.D0) CYCLE
            IF (.NOT.SET) THEN
               FKAA = SLATER(IA,IA,IA,IA,K)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)
            EMT(IDIAG) = EMT(IDIAG) + COEF*FKAA
         END DO
      END DO
      IF (K < KM) GO TO 6

!=======================================================================
!                    0
!   Piece involving F (a,b) integrals
!=======================================================================

      DO IA = 1, NW - 1
         DO IB = IA + 1, NW
            SET = .FALSE.
            DO IR = MYID + 1, NCF, NPROCS
               COEF = FCO(0,IR + NCFPAT,IA,IB)
               IF (COEF == 0.D0) CYCLE
               IF (.NOT.SET) THEN
                  F0AB = SLATER(IA,IB,IA,IB,0)
                  SET = .TRUE.
               ENDIF
            ! IDIAG = IENDC(IR-1)+1
               IDIAG = IENDC(IR)
               EMT(IDIAG) = EMT(IDIAG) + COEF*F0AB
            END DO
         END DO
      END DO

!=======================================================================
!                    k
!   Piece involving G (a,b) integrals
!=======================================================================

      KM = 0
      K = -1
   12 CONTINUE
      K = K + 1
      DO IA = 1, NW - 1
         NKJIA = NKJ(IA)
         DO IB = IA + 1, NW
            NKJIB = NKJ(IB)
            SET = .FALSE.
            IF (NAK(IA)*NAK(IB) > 0) THEN
               KMIN = ABS((NKJIA - NKJIB)/2)
            ELSE
               KMIN = ABS((NKJIA - NKJIB)/2) + 1
            ENDIF
            IF (MOD(K - KMIN,2) /= 0) CYCLE

            KMAX = (NKJIA + NKJIB)/2
            KM = MAX0(KMAX,KM)
            IF (K<KMIN .OR. K>KMAX) CYCLE

            DO IR = MYID + 1, NCF, NPROCS
               COEF = GCO(K,IR + NCFPAT,IA,IB)
               IF (COEF == 0.D0) CYCLE
               IF (.NOT.SET) THEN
                  GKAB = SLATER(IA,IB,IB,IA,K)
                  SET = .TRUE.
               ENDIF
               ! IDIAG = IENDC(IR-1)+1
               IDIAG = IENDC(IR)
               EMT(IDIAG) = EMT(IDIAG) + COEF*GKAB
            END DO
         END DO
      END DO
      IF (K < KM) GO TO 12

!=======================================================================
!   Local storage for reading mcpXXX files
!=======================================================================

      NDIM = 1
      CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETHAM' )
      CALL ALLOC (INDX, NDIM, 'INDX', 'SETHAM')

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

      NCONTR_tot2(1) = NCONTR_tot2(1) + 1
      CALL read1_mem(31,NCONTR_tot2(1),LAB, NCONTR)
      IF (IOS /= 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR'
      DO WHILE(LAB/=0 .OR. NCONTR/=0)

         !*** decode the label of I(ab)
         IA = MOD(LAB,KEY)
         IB = LAB/KEY

         !*** Compute radial integral I(ab)
         TEGRAL = RINTI(IA,IB,0)

         ! Read column index, sparse matrix index, and coefficient
         ! for all contributions from this integral.
         IF (NCONTR > NDIM) THEN
            CALL DALLOC (COEFF, 'COEFF', 'SETHAM')
            CALL DALLOC (INDX, 'INDX', 'SETHAM')
            NDIM = NCONTR
            CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETHAM')
            CALL ALLOC (INDX, NDIM, 'INDX', 'SETHAM')
         ENDIF
         CALL read2_mem(31,NCONTR_tot(1),NCONTR,ICLMNDUM)
         NCONTR_tot(1) = NCONTR_tot(1) + NCONTR

         !*** Store all the contributions from this integral
         DO I = 1, NCONTR
            LOC = INDX(I)
            LOD = LOC   !  (Convert type)
            IF (LOD > NELMNT) THEN
               WRITE (6, *) '  Error in computing 1-e contribution'
               WRITE (6, *) '  LOC = ', LOC , '  NELMNT = ', NELMNT
               STOP
            ENDIF
            EMT(LOC) = EMT(LOC) + TEGRAL*COEFF(I)
         END DO

         NCONTR_tot2(1) = NCONTR_tot2(1) + 1
         CALL read1_mem(31,NCONTR_tot2(1),LAB, NCONTR)
         IF (IOS == 0) CYCLE
         STOP 'IOS .NE. 0 when reading LAB, NCONTR'
      END DO

      DO NFILE = 32, 32 + KMAXF
         K = NFILE - 32

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

         NCONTR_tot2(NFILE-30) = NCONTR_tot2(NFILE-30) + 1
         CALL read1_mem(NFILE,NCONTR_tot2(NFILE-30),LAB, NCONTR)
         IF (IOS /= 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR 2'
         DO WHILE(LAB/=0 .OR. NCONTR/=0)

            !                         k
            !*** decode the label of R (abcd)
            ID = MOD(LAB,KEY)
            LAB = LAB/KEY
            IB = MOD(LAB,KEY)
            LAB = LAB/KEY
            IC = MOD(LAB,KEY)
            IA = LAB/KEY

            !*** Compute radial integral
            TEGRAL = SLATER(IA,IB,IC,ID,K)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            IF (NCONTR > NDIM) THEN
               CALL DALLOC (COEFF, 'COEFF', 'SETHAM')
               CALL DALLOC (INDX, 'INDX', 'SETHAM')
               NDIM = NCONTR
               CALL ALLOC (COEFF, NDIM, 'COEFF', 'SETHAM')
               CALL ALLOC (INDX, NDIM, 'INDX', 'SETHAM')
            ENDIF
            CALL read2_mem(NFILE,NCONTR_tot(NFILE-30),NCONTR,ICLMNDUM)
            NCONTR_tot(NFILE-30) = NCONTR_tot(NFILE-30) + NCONTR

            !*** Store all the contributions from this integral
            DO I = 1, NCONTR
               LOC = INDX(I)
               LOD = LOC
               IF (LOD > NELMNT) THEN   !! NELMENT (from Hmat_C) is INT*8
                  WRITE (6, *) '  Error in computing 2-e contribution'
                  WRITE (6, *) '  LOC= ', LOC , '  NELMNT = ', NELMNT
                  STOP
               ENDIF
               EMT(LOC) = EMT(LOC) + TEGRAL*COEFF(I)
            END DO

            NCONTR_tot2(NFILE-30) = NCONTR_tot2(NFILE-30) + 1
            CALL read1_mem(NFILE,NCONTR_tot2(NFILE-30),LAB, NCONTR)
            IF (IOS == 0) CYCLE
            STOP 'IOS .NE. 0 when reading LAB, NCONTR 2'
         END DO
      END DO

!=======================================================================
!   Deallocate local storage
!=======================================================================

      CALL DALLOC (COEFF, 'COEFF', 'SETHAM')
      CALL DALLOC (INDX, 'INDX', 'SETHAM')

      RETURN
      END SUBROUTINE SETHAM
