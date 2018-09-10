!***********************************************************************
!                                                                      *
      SUBROUTINE MCPOUT(NAME, IK, NTESTG, INPCI) 
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                                                      *
!                   T  (ab)                                            *
!                    rs                                                *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   are orbital sequence numbers.  r and s are configuration           *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
!               [GENMCP]: QSORT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: JUne 1996     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW, KEYORB
      USE memory_man
      USE FOPARM_C 
      USE MCP_C 
      USE PRNT_C 
      USE SYMA_C 
      USE STAT_C 
      USE BLK_C 
      USE orb_C,           ONLY: ncf, nw,nak, iqa
      USE jqjc_C
      USE orbord_C
      USE mcpdata_C
      USE sbdat_C,         ONLY: KAMAX,NLMAX,NAKINVII,NAKINVFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE onescalar_I
      USE cord_I 
      USE itjpo_I 
      USE angdata_I 
      USE qqsort_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IK 
      INTEGER, INTENT(IN) :: NTESTG 
      INTEGER  :: INPCI 
      CHARACTER  :: NAME*24 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NF = 200 
      INTEGER, PARAMETER :: NVMAX = 100 
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEYSQ = KEY*KEY 
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NLMAX) :: LLISTT 
      INTEGER, DIMENSION(NNNW) :: NAKINV 
      INTEGER, DIMENSION(NLMAX) :: NSHL, NINL 
      INTEGER :: I, NTESTL, NTEST, KA, IOPAR, J, NCF0, IBLK, K, &
         JA, JB, IBB, IA, IB, LAB, JAN, JBN, L 
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL 
      REAL(DOUBLE), DIMENSION(NVMAX) :: SC 
      REAL(DOUBLE), DIMENSION(20*NLMAX*NLMAX) :: CIROT 
      REAL(DOUBLE), DIMENSION(10000) :: EVSC 
      LOGICAL   :: F0INT, LINCR, RESTRT, COMP, AVAIL 
      CHARACTER :: CNUM*20, CK*2 
!-----------------------------------------------
!
! Locals
!     POINTER (pscr, scr(1))
!     POINTER (pciout, ciout(1))
!
 
      WRITE (6, *) 'NBLOCK,(NCFBLK(i),i=1,NBLOCK)' 
      WRITE (6, *) NBLOCK, (NCFBLK(I),I=1,NBLOCK) 
 
      NTESTL = 0 
      NTEST = MAX(NTESTG,NTESTL) 
      NTEST = 0 
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      KA = 0 
      IOPAR = 1 
!
!   Check if angular data is available. If available read this data.
!   If not available calculate the data
!
      CALL ANGDATA (NAME, AVAIL, KAMAX) 
      WRITE (6, *) 'AVAIL=', AVAIL 
      IF (AVAIL) RETURN  
 
      IF (IK == 1) THEN 
         NAKINV(:NNNW) = NAKINVII(:NNNW) 
      ELSE 
         NAKINV(:NNNW) = NAKINVFF(:NNNW) 
      ENDIF 
 
 
 
      WRITE (6, *) ' open sorted ang. file .TB(NF)', NF 
      J = INDEX(NAME,' ') 
      OPEN(UNIT=NF, FILE=NAME(1:J-1)//'.TB', STATUS='UNKNOWN', FORM=&
         'UNFORMATTED', POSITION='asis') 
 
      NCF0 = 1 
      DO IBLK = 1, NBLOCK 
!
!   Open scratchfiles to dump the T coefficients for each kappa
!
         DO K = 1, KAMAX 
            OPEN(UNIT=80 + K, STATUS='UNKNOWN', FORM='UNFORMATTED', POSITION=&
               'asis') 
         END DO 
!
!   Initialize the counters for the total number of T coefficients
!
         LLISTT(:NLMAX) = 0 
!
!   JA and JB respectively refer to the initial and final states
!   in the list of NCF configurations
!
         DO JA = NCF0, NCFBLK(IBLK) 
            IF (MOD(JA,100)==0 .AND. IK==1) THEN 
               WRITE (*, *) ' JA1 =', JA, JA - NCF0 + 1 
            ELSE IF (MOD(JA,100)==0 .AND. IK==2) THEN 
               WRITE (*, *) ' JA2 =', JA, JA - NCF0 + 1 
            ENDIF 
!
            DO JB = NCF0, NCFBLK(IBLK) 
!
!   Call the MCT package to compute T coefficients
!
               IF (NTRANS == 1) THEN 
                  COMP = .FALSE. 
                  IF (IK == 1) THEN 
                     DO IBB = 1, JQJ1 
                        IF (ITJPO(JA) /= ITJQJ1(IBB)) CYCLE  
                        COMP = .TRUE. 
                     END DO 
                  ELSE 
                     DO IBB = 1, JQJ2 
                        IF (ITJPO(JA) /= ITJQJ2(IBB)) CYCLE  
                        COMP = .TRUE. 
                     END DO 
                  ENDIF 
               ELSE 
                  COMP = .TRUE. 
               ENDIF 
 
 
               IF (.NOT.COMP) CYCLE  
!            write(*,*) JA,JB
!GG               CALL TNSRJJ (KA, IOPAR, JA, JB, IA, IB, TSHELL) 
              CALL ONESCALAR(JA,JB,IA,IB,TSHELL)
               IF (IA == 0) CYCLE  
               IF (IA == IB) THEN 
                  DO IA = 1, NW 
!
!   If T coefficient is greater than zero and the kappa quantum numbers
!   of the two orbitals are the same dump to file
!   In a later version use a buffer with a reasonable record length
!
                     IF (DABS(TSHELL(IA)) <= CUTOFF) CYCLE  
                     LLISTT(NAKINV(IA)) = LLISTT(NAKINV(IA)) + 1 
                     LAB = IA*KEY + IA 
                     WRITE (80 + NAKINV(IA)) JA - NCF0 + 1, JB - NCF0 + 1, LAB&
                        , TSHELL(IA) 
                  END DO 
               ELSE 
                  IF (DABS(TSHELL(1))>CUTOFF .AND. NAK(IA)==NAK(IB)) THEN 
                     LLISTT(NAKINV(IA)) = LLISTT(NAKINV(IA)) + 1 
                     IF (NORDII==0 .AND. NORDFF==0) THEN 
!
!   Experssion for normal orbital ordering
!
                        LAB = IA*KEY + IB 
                        JAN = JA - NCF0 + 1 
                        JBN = JB - NCF0 + 1 
                     ELSE IF (NORDII==1 .AND. NORDFF==1) THEN 
!
!   Experssion for reversed orbital ordering
!
                        LAB = IB*KEY + IA 
                        JAN = JB - NCF0 + 1 
                        JBN = JA - NCF0 + 1 
                     ELSE 
                        WRITE (*, *) 'SOMETHING WRONG'
                        STOP
                     ENDIF 
                     WRITE (80 + NAKINV(IA)) JAN, JBN, LAB, TSHELL(1) 
                  ENDIF 
               ENDIF 
!
            END DO 
         END DO 
 
 
!
!   sort the MCP data into inegral based lists.
!
         DO L = 1, KAMAX 
            IF (LLISTT(L) > 0) THEN 
               CALL QQSORT (L, LLISTT(L), IK, NAME, KAMAX) 
            ELSE 
               IF (L == 1) WRITE (NF) NCF, NW, KAMAX 
               WRITE (NF) LLISTT(L), LLISTT(L) 
            ENDIF 
         END DO 
!
!  Close the angular files
!
         DO L = 1, KAMAX 
            CLOSE(L + 80, STATUS='DELETE') 
         END DO 
!
         NCF0 = NCFBLK(IBLK) + 1 
      END DO 
!
!   Deallocate storage that is no longer required. This was
!   allocated in lodcsl.
!
      CALL DALLOC (IQA, 'IQA', 'MCPOUT') 
      CALL DALLOC (JQSA, 'JQSA', 'MCPOUT') 
      CALL DALLOC (JCUPA, 'JCUPA', 'MCPOUT') 
      RETURN  
      END SUBROUTINE MCPOUT 
