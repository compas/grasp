!***********************************************************************
!                                                                      *
      SUBROUTINE SETLAGmpi(EOL)
!                                                                      *
!   Sets up the data structure  pertaining to the Lagrange multipli-   *
!   ers  on the first entry;  on subsequent calls it  determines new   *
!   estimates for the multipliers.                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
!               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
!   MPI version by Xinghong He              Last update: 03 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:25:08   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNP
      USE memory_man
      USE ORBA_C
      USE core_C
      USE def_C
      USE fixd_C
      USE grid_C, ONLY: n, rpor
      USE lagr_C
      USE orb_C
      USE mpi_C
      USE pote_C, ONLY: yp, xp, xq
      USE scf_C
      USE tatb_C, ONLY: mtp, ta
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcof_I
      USE ypot_I
      USE xpot_I
      USE dacon_I
      USE quad_I
      USE rinti_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL  :: EOL
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: P001 = 1.0D-01
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITWICE, LIRAW, LI, LIP1, NAKLJ, LJRAW, LJ, IECCLI, L1, L2, &
         JLAST, MLAST, M, J, I
      REAL(DOUBLE), DIMENSION(NNNP) :: YPJ, YPM, XPJ, XPM, XQJ, XQM
      REAL(DOUBLE) :: EPS,UCFJ,UCFM,RESULT,RIJM,QDIF,OBQDIF,OBQSUM,TMP
!cff  March 2019
      LOGICAL :: FIRST, FIXLI, FIXLJ, FULLI, FULLJ, VLI, VLJ
!-----------------------------------------------
!
      DATA FIRST/ .TRUE./
!
!-----------------------------------------------------------------------

      IF (FIRST) THEN

!=======================================================================
!   Determine the total number of Lagrange multipliers and store
!   their indeces in IECC(1:NEC). Memories are allocated for IECC
!   and ECV. The only outputs are NEC and IECC(). This is done only
!   for the first time this routine is called.
!
!   Implementation changes:
!   The "DO LIraw ..." is executed twice (if NEC obtained in the
!   first time > 0) so the memory allocation is done only once.
!
!   This part is not distributed.
!=======================================================================

         EPS = ACCY*0.01D0      ! criterion to see if an orb is occupied
         DO ITWICE = 1, 2
            NEC = 0
         DO LJraw = NCORE+1, NW
               LJ = iorder(LJraw)
            NAKLJ = NAK(LJ)
            VLJ = .NOT. LFIX(LJ)
               FULLJ = ABS ( UCF(LJ)-DBLE (NKJ(LJ)+1) ) .LT. EPS
            DO LIraw = 1, LJraw-1
               LI = iorder(LIraw)
               VLI = .NOT. LFIX(LI)
               FULLI = ABS ( UCF(LI)-DBLE (NKJ(LI)+1) ) .LT. EPS
               IF  (NAK(LI) .EQ. NAKLJ) then
                  If  (VLI  .OR. VLJ ) then          !at least one varid
!                                                    ! but not (both varied and full)
                     If (.NOT. ((VLI .AND. VLJ) .AND. (FULLI .AND. FULLJ))) THEN 
                  NEC = NEC + 1
                  !*** Encode index at 2nd round ***
                  IF (itwice == 2) IECC(NEC) = LI + KEY * LJ
               ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO


            IF (ITWICE==1 .AND. NEC>0) THEN
               CALL ALLOC (ECV, NEC, 'ECV', 'SETLAGmpi')
               CALL ALLOC (IECC, NEC, 'IECC', 'SETLAGmpi')
            ELSE
               EXIT
            ENDIF
         END DO                                  !itwice

!=======================================================================
!   Print information about Lagrange multipliers
!=======================================================================

         IF (MYID == 0) THEN
            IF (NEC == 0) THEN
               WRITE (*, 302)
            ELSE
               WRITE (*, 304)
               DO LI = 1, NEC
               !*** Decode index ***
                  IECCLI = IECC(LI)
                  L1 = IECCLI/KEY
                  L2 = IECCLI - KEY*L1
                  WRITE (*, 305) NP(L2), NH(L2), NP(L1), NH(L1)
               END DO
            ENDIF
         ENDIF
         FIRST = .FALSE.
      ENDIF

!FF+GG  12/07/05
!     Lagrange multipliers need to be computed also on the first call
!     RETURN

      IF (NEC == 0) RETURN
      IF (MYID == 0) WRITE (*, 306)
      JLAST = 0
      MLAST = 0

      DO LI = 1, NEC
         !*** Decode index ***
         IECCLI = IECC(LI)
         M = IECCLI/KEY
         J = IECCLI - KEY*M
!
         IF (J /= JLAST) THEN
            UCFJ = UCF(J)
            CALL SETCOF (EOL, J)
            CALL YPOT (J)
            CALL XPOT (J)
            CALL DACON
            YPJ(:N) = YP(:N)
            XPJ(:N) = XP(:N)
            XQJ(:N) = XQ(:N)
            JLAST = J
         ENDIF
!
         IF (M /= MLAST) THEN
            UCFM = UCF(M)
            CALL SETCOF (EOL, M)
            CALL YPOT (M)
            CALL XPOT (M)
            CALL DACON
            YPM(:N) = YP(:N)
            XPM(:N) = XP(:N)
            XQM(:N) = XQ(:N)
            MLAST = M
         ENDIF
!
         MTP = MAX(MF(J),MF(M))
!
         IF (LFIX(M)) THEN
            TA(1) = 0.D0
            DO I = 2, MTP
               TA(I) = RPOR(I)*((PF(I,M)*XQJ(I)-QF(I,M)*XPJ(I))*C+       &
                                (PF(I,M)*PF(I,J)+QF(I,M)*QF(I,J))*YPJ(I))
            END DO

            CALL QUAD (RESULT)
            RIJM = RINTI(M,J,1)
            ECV(LI) = (RESULT - RIJM / nprocs)*UCFJ
! start dbg
!           WRITE (81,*)'1, RESULT, RIJM, UCFJ, ECV, TA' ! dbg
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV ! dbg
!           DO i = 1, MTP ! dbg
!              WRITE (81,*) i, TA(i), r(i), rp(i) ! dbg
!           ENDDO ! dbg
! end dbg

         ELSE IF (LFIX(J)) THEN
            TA(1) = 0.D0
            DO I = 2, MTP
               TA(I) = RPOR(I)*((PF(I,J)*XQM(I)-QF(I,J)*XPM(I))*C+      &
                                (PF(I,J)*PF(I,M)+QF(I,J)*QF(I,M))*YPM(I))
            END DO

!start dbg
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i)
!               WRITE (83,*) i, r(i), rp(i), rpor(i)
!               WRITE (84,*) i, pf(i,j), qf(i,j)
!               write(85,*) i,ypm(i)
!              write(86,*)i,xpm(i),xqm(i)
!           ENDDO
! end dbg
            CALL QUAD (RESULT)

            RIJM = RINTI(J,M,1)                  !/ nprocs
            ECV(LI) = (RESULT - RIJM / nprocs)*UCFM
!start dbg
!           WRITE (81,*)'2, RESULT, RIJM, UCFM, ECV, TA'
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV, r(i), rp(i)
!end dbg


         ELSE
               OBQSUM = 1.D0/UCFJ + 1.D0/UCFM
               TA(1) = 0.D0
               DO I = 2, MTP
                  TA(I) = RPOR(I)*((PF(I,M)*XQJ(I)-QF(I,M)*XPJ(I)         &
                                   +PF(I,J)*XQM(I)-QF(I,J)*XPM(I))*C      &
                         +(YPJ(I)+YPM(I))*(PF(I,M)*PF(I,J)+QF(I,M)*QF(I,J)))
               END DO

               CALL QUAD (RESULT)
               RIJM = RINTI(M,J,1)               !/ nprocs
               ECV(LI) = (RESULT - 2.D0*RIJM / nprocs)/OBQSUM
!start dbg
!           WRITE (81,*)'4, RESULT, RIUJM, OBQSUM, ECV, TA'
!           WRITE (81,*)RESULT, RIUJM, OBQSUM, ECV
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i), r(i), rp(i)
!           ENDDO
!end dbg

!           ENDIF 
         ENDIF
!=======================================================================
!  Collect contributions from all nodes.
!  Another alternative is to modify mcpmpi to let every node
!  have the same set of "LAB" and "NCONTR". So an easier aggregation
!  can be done for arrays DA, XA, YA in setcof. And so the computation
!  of YPOT, XPOT, DACON can be distributed.
!=======================================================================
         CALL MPI_Allreduce (ECV(LI),tmp,1,MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
         ECV(LI) = tmp

         IF(MYID == 0)WRITE (*, 307) NP(J), NH(J), NP(M), NH(M), ECV(LI)

      END DO

!db      close(81)
!db      close(82)


  302 FORMAT(/,'Lagrange multipliers are not required')
  304 FORMAT(/,'Include Lagrange multipliers between:'/)
  305 FORMAT(13X,2(2X,1I2,1A2))
  306 FORMAT(/,'Lagrange multipliers:'/)
  307 FORMAT(13X,2(2X,1I2,1A2),2X,1P,D16.9)

      RETURN
      END SUBROUTINE SETLAGmpi
