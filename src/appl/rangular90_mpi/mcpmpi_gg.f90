!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      SUBROUTINE mcpmpi (nb, RESTRT, fhead)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                       k                              *
!                   T  (ab)            V  (abcd)                       *
!                    rs                 rs                             *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   c and d are orbital sequence numbers.  r and s are configuration   *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, RKCO_GG,       *
!                        TNSRJJ.                                       *
!               [GENMCP]: FNDBEG, SETSDA, SORT.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
!   Modified by C. Froese Fischer for block computation.               *
!   Modified by G. Gaigalas and J. Bieron for new spin-angular         *
!   integration                                        01 April 2012   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW, KEYORB
      USE memory_man
      USE mpi_C,           ONLY:  myid, nprocs
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outsdampi_I
      USE fndbeg_I
      USE alcbuf_I
      USE rkco_GG_I
      USE setsda_I
      USE sort_I
!-----------------------------------------------
!   C O M M O N  B L O C K S
!-----------------------------------------------
      USE BUFFER_C,  ONLY: NVCOEF, LABEL, COEFF
      USE DEBUG_C,   ONLY: LDBPA
      USE DEFAULT_C, ONLY: NDEF
      USE iccu_C,    ONLY: ICCUT
      USE MCP_C,     ONLY: KMAX, DIAG, LFORDR
      USE ORB_C,     ONLY: NCF, IQA
      USE STAT_C,    ONLY: JQSA, JCUPA
      IMPLICIT NONE
      EXTERNAL CORD
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NB
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      LOGICAL , INTENT(IN) :: RESTRT
      CHARACTER(len=*), INTENT(IN) :: FHEAD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!GG      REAL(DOUBLE), PARAMETER :: CUTOFF = 0.0D0
!cjb  REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-08
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(:), pointer :: LLISTV
      INTEGER :: JASTRT,JBSTRT,NPOS,LLISTT,K,JA,JB,IA,IB,LAB,I,IC, &
         ID, NSWAP, ISWAP, LAC, LBD, NTGI
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF
      LOGICAL :: F0INT, LINCR
!-----------------------------------------------
!
!cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD '
!cjb           PRINT *, ' mcpmpi_gg.f90 FHEAD = ', FHEAD
!
!   Allocate storage to the array that stores the lengths
!   of the lists of V coefficients of each multipolarity
!
      CALL ALLOC(LLISTV, 0, kmax, 'LLISTV', 'MCPMPI')
!
!   If this is a restart, determine the pair of CSFs with which
!   the computation should begin.
!   If not restart, initialize the counters for the total number of
!   T coefficients and V coefficients of each multipolarity.
!
      IF (RESTRT) THEN
         PRINT *, 'Not ready for RESTART'
         STOP
         CALL FNDBEG (JASTRT, JBSTRT, npos, LLISTT, LLISTV)
      ELSE
         JASTRT = 1 + myid
         npos = 0
         LLISTT = 0
         DO K = 0, KMAX
            LLISTV(K) = 0
         ENDDO
      ENDIF
!
!   Allocate storage for the arrays in BUFFER
!
      CALL ALCBUF (1)
!
!   Write header to  .dbg  file if appropriate
!
      IF (LDBPA(2) .OR. LDBPA(3)) WRITE (99,300)
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
      DO 5 JA = JASTRT, NCF, nprocs
         JBSTRT = 1
         DO 4 JB = JBSTRT, JA
            IF (DIAG .OR. (LFORDR .AND. (JB .GT. ICCUT(nb))) ) THEN
               IF (JB /= JA) CYCLE
            END IF
!
!   LINCR is .TRUE. if npos is to be incremented by 1; there
!   is always a diagonal element in each column
!
            IF (JB /= JA) THEN
               LINCR = .TRUE.
            ELSE
               npos = npos + 1
               LINCR = .FALSE.
            ENDIF
            IF (JB /= JA) THEN
!   Compute T coefficients
               CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
               IF (IA /= 0 .AND. IA /= IB .AND.       &
                      ABS (TSHELL(1)) > CUTOFF) THEN
                  npos = npos + 1
                  LINCR = .FALSE.
                  LLISTT = LLISTT + 1
                  LAB = MIN (IA,IB) * KEY + MAX (IA,IB)
                  WRITE (31) JA, npos, LAB, TSHELL(1)
               ENDIF
            ENDIF
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
            NVCOEF = 0
!            CALL RKCO_GG (JA, JB, 0, 1)
            CALL RKCO_GG (JA, JB, CORD, 0, 1)
            DO 2 I = 1, NVCOEF
               VCOEFF = COEFF(I)
               IF (ABS (VCOEFF) > CUTOFF) THEN
                  IA = LABEL(1,I)
                  IB = LABEL(2,I)
                  IC = LABEL(3,I)
                  ID = LABEL(4,I)
                  K  = LABEL(5,I)
                  F0INT=(K.EQ.0).AND.(IA.EQ.IC).AND.(IB.EQ.ID)
                  IF (.NOT. F0INT) THEN
! Swap index to make sure IA <= IC, IB <= ID and record the number
! of swaps
                     NSWAP = 0
                     IF (IA > IC) THEN
                        ISWAP = IC
                        IC = IA
                        IA = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (IB > ID) THEN
                        ISWAP = ID
                        ID = IB
                        IB = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (LINCR) THEN
                        npos = npos + 1
                        LINCR = .FALSE.
                     ENDIF
                     LLISTV(K) = LLISTV(K) + 1
                     LAC = IA * KEY + IC
                     LBD = IB * KEY + ID
                     IF (LAC .LT. LBD) THEN
                        LAB = LAC * KEYSQ + LBD
                     ELSE
                        LAB = LBD * KEYSQ + LAC
                     ENDIF
                     WRITE (32+K) JA, npos, LAB, VCOEFF
                  ENDIF
               ENDIF
    2       CONTINUE
!
!   All angular coefficients for this pair of CSFs have been
!   generated; update file 30
!
            IF ((JB .EQ. JA) .OR. (.NOT. LINCR))   &
               WRITE (30) JA, JB, npos
    4    CONTINUE
         if (mod(JA-1,100) == 0) PRINT '(A4,I7,A9,I10,A9,I3,A8,I3)',  &
            'Row ',JA,' nnonz = ',npos,' block = ',nb,' myid = ',myid
    5 CONTINUE
!
!   Deallocate storage that is no longer required
!
      CALL DALLOC (iqa, 'IQA', 'MCPMPI')
      CALL DALLOC (jqsa, 'JQSA', 'MCPMPI')
      CALL DALLOC (jcupa, 'JCUPA', 'MCPMPI')
      CALL ALCBUF (3)
!
!   Write out a report for this run
!
      IF (NDEF .NE. 0 .AND. myid .EQ. 0) THEN
         WRITE (24,*)
         WRITE (24,*) LLISTT, ' T coefficients generated;'
         DO K = 0, KMAX
            WRITE (24,*) LLISTV(K),' V(k=',K,') coefficients generated;'
         ENDDO
      ENDIF
!
!   Set up sparse structure definition arrays in file 30
!

      REWIND (30)
      CALL SETSDA (outsdampi, npos, LDBPA(4), nb, myid, nprocs, fhead)
      CLOSE (30, STATUS = 'DELETE')
!
!   Sort MCP coefficients into integral-based lists
!
!   1) the T coefficient
      REWIND (31)
      CALL SORT (31, (LLISTT), ntgi, LDBPA(2), nb, fhead)
      IF (NDEF .NE. 0 .AND. myid .EQ. 0)      &
         WRITE (24,*) ntgi, ' I(ab) integrals;'
      CLOSE (31, STATUS = 'DELETE')
!   2) the V coefficient
      DO k = 0, KMAX
         REWIND (k+32)
         CALL SORT (k+32, (LLISTV(k)), ntgi, LDBPA(3), nb, fhead)
         IF (NDEF .NE. 0 .AND. myid .EQ. 0)       &
            WRITE (24,*) ' k = ', k, ': ', ntgi, ' Slater integrals;'
         CLOSE (k+32, STATUS = 'DELETE')
      ENDDO
       CALL DALLOC(LLISTV, 'LLISTV', 'MCPMPI')
  300 FORMAT (/'From MCPMPI:')
      RETURN
      END SUBROUTINE mcpmpi
