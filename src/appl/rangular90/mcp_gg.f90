!***********************************************************************
!                                                                      *
      SUBROUTINE mcp (nb, RESTRT, myid, nprocs, fhead)
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
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outsda_I
      USE fndbeg_I
      USE alcbuf_I
      USE rkco_GG_I
      USE setsda_I
      USE sort_I
      USE sortmem_I
      USE allocCheck_I
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
      USE sacoef_C
      IMPLICIT NONE
      EXTERNAL CORD
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NB
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      LOGICAL, INTENT(IN) :: RESTRT
      CHARACTER(len=*), INTENT(IN) :: FHEAD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEYORB*KEYORB
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! GG begin
!     2^(23)
      INTEGER, PARAMETER :: NGGMAX = 8388608
!     2^(24)
!        PARAMETER (NGGMAX = 16777216)
!     2^(27)
!        PARAMETER (NGGMAX   = 134217728)
!        PARAMETER (NGGMAX = 536870912)
! GG end
      INTEGER, DIMENSION(:), pointer :: LLISTV, MARK
      INTEGER :: JASTRT,JBSTRT,NPOS,LLISTT,K,JA,JB,IA,IB,LAB,I,IC, &
         ID, NSWAP, ISWAP, LAC, LBD, NTGI
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE) :: VCOEFF
      LOGICAL :: F0INT, LINCR
      INTEGER :: IREZ, IGGMAX, IGGMAX_K, MARKT, NEWSIZ, JJ
!-----------------------------------------------
!
!   Allocate storage to the array that stores the lengths
!   of the lists of V coefficients of each multipolarity
!
! GG begin
      IREZ = 1
      IGGMAX = NGGMAX
      IGGMAX_K = IGGMAX/(2+KMAX)
      IGGMAX = IGGMAX_K*(2+KMAX)
!GG      PRINT*, "The program alloc Memory for sorting:",IGGMAX
      CALL ALLOC(LLISTV, 0, kmax, 'LLISTV', 'MCP')
      CALL ALLOC (MARK, KMAX + 1,'MARK', 'MCP')
      CALL ALLOC (TCOEFF,IGGMAX_K,KMAX+2,'TCOEFF','MCP')
      CALL ALLOC (IICLMN,IGGMAX_K,KMAX+2,'TICLMN','MCP')
      CALL ALLOC (IINDEX,IGGMAX_K,KMAX+2,'TINDEX','MCP')
      CALL ALLOC (ILABEL,IGGMAX_K,KMAX+2,'TLABEL','MCP')
      MARKT = 0
      MARK = 0
! GG end
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
               IF (JB.NE.JA) CYCLE
            END IF
!   LINCR is .TRUE. if npos is to be incremented by 1; there
!   is always a diagonal element in each column
            IF (JB .NE. JA) THEN
               LINCR = .TRUE.
            ELSE
               npos = npos + 1
               LINCR = .FALSE.
            ENDIF
            IF (JB .NE. JA) THEN
               CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
               IF (IA .NE. 0 .AND. IA .NE. IB .AND.       &
                      ABS (TSHELL(1)) .GT. CUTOFF) THEN
                  npos = npos + 1
                  LINCR = .FALSE.
                  LLISTT = LLISTT + 1
                  LAB = MIN (IA,IB) * KEY + MAX (IA,IB)
! GG begin
                  IF(MARKT .EQ. 0) THEN
                    IF(LLISTT .GT. IGGMAX_K) THEN
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      NEWSIZ = (IGGMAX_K+IGGMAX_K/2)*(2+KMAX)
                      CALL ALLOCCHECK (NEWSIZ,IREZ)
                      IF(IREZ .EQ. 1) THEN
                        IGGMAX = NEWSIZ
                        IGGMAX_K = IGGMAX/(2+KMAX)
                        CALL RALLOC                                    &
                             (TCOEFF,IGGMAX_K,KMAX+2,'TCOEFF','MCP')
                        CALL RALLOC                                    &
                             (IICLMN,IGGMAX_K,KMAX+2,'TICLMN','MCP')
                        CALL RALLOC                                    &
                             (IINDEX,IGGMAX_K,KMAX+2,'TINDEX','MCP')
                        CALL RALLOC                                    &
                             (ILABEL,IGGMAX_K,KMAX+2,'TLABEL','MCP')
                        IICLMN(LLISTT,1) = JA
                        IINDEX(LLISTT,1) = npos
                        ILABEL(LLISTT,1) = LAB
                        TCOEFF(LLISTT,1) = TSHELL(1)
                      ELSE
!GG                        print*,"LLISTT =",LLISTT
!GG                        print*, 
!GG     :             "The program switches to the disk version of
!sorting"
                        CALL SETTMPGG (nb, 31, 'tmp')
                        MARKT = 1
                        DO JJ = 1, LLISTT-1
                          WRITE (31) IICLMN(JJ,1), IINDEX(JJ,1),       &
                                     ILABEL(JJ,1), TCOEFF(JJ,1)
                        END DO
                        WRITE (31) JA, npos, LAB, TSHELL(1)
                      END IF
                    ELSE
                      IICLMN(LLISTT,1) = JA
                      IINDEX(LLISTT,1) = npos
                      ILABEL(LLISTT,1) = LAB
                      TCOEFF(LLISTT,1) = TSHELL(1)
                    END IF
                  ELSE
                    WRITE (31) JA, npos, LAB, TSHELL(1)
                  END IF
! GG end
               ENDIF
            ENDIF
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs. NVCOEF is initialized here but changed
!   in /rkco/cor[d] via COMMON.
!
            NVCOEF = 0
            CALL RKCO_GG (JA, JB, CORD, 0, 1)
            DO 2 I = 1, NVCOEF
               VCOEFF = COEFF(I)
               IF (ABS (VCOEFF) .GT. CUTOFF) THEN
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
                     IF (IA .GT. IC) THEN
                        ISWAP = IC
                        IC = IA
                        IA = ISWAP
                        NSWAP = NSWAP + 1
                     ENDIF
                     IF (IB .GT. ID) THEN
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
! GG begin
                  IF(MARK(K+1) .EQ. 0) THEN
                    IF(LLISTV(K) .GT. IGGMAX_K) THEN
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
                      NEWSIZ = (IGGMAX_K+IGGMAX_K/2)*(2+KMAX)
                      CALL ALLOCCHECK (NEWSIZ,IREZ)
                      IF(IREZ .EQ. 1) THEN
                        IGGMAX = NEWSIZ
                        IGGMAX_K = IGGMAX/(2+KMAX)
                        CALL RALLOC                                    &
                             (TCOEFF,IGGMAX_K,KMAX+2,'TCOEFF','MCP')
                        CALL RALLOC                                    &
                             (IICLMN,IGGMAX_K,KMAX+2,'TICLMN','MCP')
                        CALL RALLOC                                    &
                             (IINDEX,IGGMAX_K,KMAX+2,'TINDEX','MCP')
                        CALL RALLOC                                    &
                             (ILABEL,IGGMAX_K,KMAX+2,'TLABEL','MCP')
                        IICLMN(LLISTV(K),K+2) = JA
                        IINDEX(LLISTV(K),K+2) = npos
                        ILABEL(LLISTV(K),K+2) = LAB
                        TCOEFF(LLISTV(K),K+2) = VCOEFF
                      ELSE
!GG                        print*,"K=",K," LLISTV =",LLISTV(K)
!GG                        print*, 
!GG     :             "The program switches to the disk version of
!sorting"
                        CALL SETTMPGG (nb, 32+K, 'tmp')
                        MARK(K+1) = 1
                        DO JJ = 1,IGGMAX_K
                          WRITE (32+K) IICLMN(JJ,K+2),IINDEX(JJ,K+2),  &
                                       ILABEL(JJ,K+2),TCOEFF(JJ,K+2)
                        END DO
                        WRITE (32+K) JA, npos, LAB, VCOEFF
                      END IF
                    ELSE
                      IICLMN(LLISTV(K),K+2) = JA
                      IINDEX(LLISTV(K),K+2) = npos
                      ILABEL(LLISTV(K),K+2) = LAB
                      TCOEFF(LLISTV(K),K+2) = VCOEFF
                    END IF
                  ELSE
                    WRITE (32+K) JA, npos, LAB, VCOEFF
                  END IF
! GG end


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
      if (mod(ja,100) == 0) then
         IF (JA .EQ. 1 .OR. JA .EQ. NCF .OR. MOD (JA,100) .EQ. 0) THEN
            PRINT *, 'Row ', JA, '; nnonz = ', npos,';  block = ', nb
         ENDIF
      end if
    5 CONTINUE
!
!   Deallocate storage that is no longer required
!
      CALL DALLOC (iqa, 'IQA', 'MCP')
      CALL DALLOC (jqsa, 'JQSA', 'MCP')
      CALL DALLOC (jcupa, 'JCUPA', 'MCP')
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
      CALL SETSDA (outsda, npos, LDBPA(4), nb, myid, nprocs, fhead)
      CLOSE (30, STATUS = 'DELETE')
!
!   Sort MCP coefficients into integral-based lists
!
! 1) the T coefficient
      IF(MARKT .EQ. 0) THEN
         CALL SORTMEM(31,IGGMAX_K,(LLISTT), ntgi, LDBPA(2), nb, fhead)
      ELSE
         REWIND (31)
         CALL SORT (31, (LLISTT), ntgi, LDBPA(2), nb, fhead)
         CLOSE (31, STATUS = 'DELETE')
      END IF
      IF (NDEF .NE. 0 .AND. myid .EQ. 0)                   &
                       WRITE (24,*) ntgi, ' I(ab) integrals;'
! 2) the V coefficient
      DO k = 0, KMAX
        IF(MARK(K+1) .EQ. 0) THEN
         CALL SORTMEM(k+32,IGGMAX_K,(LLISTV(k)),ntgi,LDBPA(3),nb,fhead)
        ELSE
          REWIND (k+32)
          CALL SORT (k+32, (LLISTV(k)), ntgi, LDBPA(3), nb, fhead)
          CLOSE (k+32, STATUS = 'DELETE')
        END IF
          IF (NDEF .NE. 0 .AND. myid .EQ. 0)               &
            WRITE (24,*) ' k = ', k, ': ', ntgi, ' Slater integrals;'
      ENDDO
      CALL DALLOC (MARK, 'MARK','MCP')
      CALL DALLOC (TCOEFF,'TCOEFF','MCP')
      CALL DALLOC (IICLMN,'TICLMN','MCP')
      CALL DALLOC (IINDEX,'TINDEX','MCP')
      CALL DALLOC (ILABEL,'TLABLE','MCP')
      CALL DALLOC(LLISTV, 'LLISTV', 'MCP')
  300 FORMAT (/'From MCP:')
      RETURN
      END SUBROUTINE mcp
