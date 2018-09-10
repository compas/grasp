!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSLM(NCORE) 
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Updated by Xinghong He                               31 Oct 1997
!       To accept both block and non-block formats
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW
      USE memory_man
      USE debug_C
      USE def_C
      USE orb_C
      USE STAT_C 
      USE TERMS_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE prsrsl_I 
      USE convrt_I 
      USE prsrcn_I 
      USE parsjl_I 
      USE pack_I 
      USE iq_I 
      USE jqs_I 
      USE jcup_I 
      USE itjpo_I 
      USE ispar_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NW2 = 2*NNNW 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(NNNW) :: IOCC 
      INTEGER , DIMENSION(NW2) :: IQSUB 
      INTEGER , DIMENSION(NNNW) :: JX 
      INTEGER :: NCORP1, NPEEL, NPEEL2, J, NPJ, NAKJ, I, LENTH, NCFD, NREC, IOS&
         , IERR, LOC, NQS, NEWSIZ, ISPARC, NJX, IOC, IPTY, NQSN, NJXN, NPEELN, &
         NOPEN, JLAST, ILAST, IOCCI, NKJI, IFULLI, NU, JSUB, IQT, NBEG, NEND, &
         JXN, JPI, II, ITEMP, NCOREL 
      LOGICAL :: EMPTY, FULL 
      CHARACTER :: RECORD*256, RECL 
!-----------------------------------------------
!
!
!   Entry message
!
      WRITE (6, *) 'Loading Configuration Symmetry List File ...' 
!
!   Get the list of subshells
!
      NW = 0 
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      CALL PRSRSL (21, 1) 
      NCORE = NW 
      NCORP1 = NW + 1 
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      READ (21, *) 
      CALL PRSRSL (21, 2) 
      NPEEL = NW - NCORE 
      NPEEL2 = NPEEL*2 
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      DO J = NCORE + 1, NW 
         NPJ = NP(J) 
         NAKJ = NAK(J) 
         DO I = 1, NCORE 
            IF (NP(I)/=NPJ .OR. NAK(I)/=NAKJ) CYCLE  
            WRITE (6, *) 'LODCSL: The lists of core and' 
            WRITE (6, *) ' peel subshells must form' 
            WRITE (6, *) ' disjoint sets.' 
            STOP  
         END DO 
      END DO 
!
!   Print the number of relativistic subshells
!
      IF (NW > 1) THEN 
         CALL CONVRT (NW, RECORD, LENTH) 
         WRITE (6, *) ' there are '//RECORD(1:LENTH)//&
            ' relativistic subshells;' 
      ELSE 
         WRITE (6, *) ' there is 1 relativistic subshell;' 
      ENDIF 
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      CALL DALLOC (IQA,   'IQA',   'LODCSLM') 
      CALL DALLOC (JQSA,  'JQSA',  'LODCSLM') 
      CALL DALLOC (JCUPA, 'JCUPA', 'LODCSLM') 
      NCFD = 2 
      CALL ALLOC (IQA,   NNNW,    NCFD, 'IQA',   'LODCSLM') 
      CALL ALLOC (JQSA,  NNNW, 3, NCFD, 'JQSA',  'LODCSLM') 
      CALL ALLOC (JCUPA, NNNW,    NCFD, 'JCUPA', 'LODCSLM') 
!
!   Skip the header for the list of CSFs
!
      READ (21, *) 
!
!   NREC is the sequence number of the last record read in the
!   Configuration Symmetry List File
!
      NREC = 5 
!
!   There must be three records for each CSF: For instance,
!
!    4s ( 2) 4p-( 2) 4p ( 3) 4d-( 2) 4d ( 5) 4f-( 6) 4f ( 4)
!                        3/2       0     5/2             2;4
!                                           1               3-
!   Each record is as follows:
!      (1) Peel subshell occupation number (q) specification.
!      (2) Subshell total angular momentum quantum number
!          (v_, J_sub) specifications. Subshell total angular
!          momentum quantum numbers (J_sub) must be specified
!          for all open subshells, even if this can be
!          deduced from the subshell occupation. Seniority
!          quantum numbers (v) must be specified for
!          subshells with j = 7/2, q = 4, J_sub = 2 or 4.
!          The seniority quantum number must precede the
!          total angular momentum quantum number: in the
!          example above, for the 4f_7/2 subshell, q = 4
!          and J_sub = 4, whence it is necessary to specify
!          v --- 2 in this case.
!      (3) The minimum number of intermediate angular
!          momentum quantum numbers (X) as well as the final
!          angular momentum quantum number (J) immediately
!          followed by the sign of the parity (P) must be
!          specified on this record. In the example above,
!          the J_sub = 0 for the 4d_3/2 subshell, whence
!          it is unnecessary to specify its coupling to all
!          preceding subshells.
!   These conventions have been chosen so as to render the CSF
!   specifications easily interpreted by the user
!
      NCF = 0 
    3 CONTINUE 
      NCF = NCF + 1 
!
      READ (21, '(A)', IOSTAT=IOS) RECORD 
!**********************************************************************
!blk*
!   To skip the border line added to mark the end of a block
!
      IF (RECORD(1:2) == ' *') READ (21, '(A)', IOSTAT=IOS) RECORD 
!**********************************************************************
 
      IF (IOS == 0) THEN 
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         CALL PRSRCN (RECORD, NCORE, IOCC, IERR) 
         IF (IERR /= 0) GO TO 26 
!
!   Read the J_sub and v quantum numbers
!
         READ (21, '(A)', IOSTAT=IOS) RECORD 
         IF (IOS /= 0) THEN 
            WRITE (6, *) 'LODCSL: Expecting subshell quantum' 
            WRITE (6, *) ' number specification;' 
            GO TO 26 
         ENDIF 
         LOC = LEN_TRIM(RECORD) 
         CALL PARSJL (1, NCORE, RECORD, LOC, IQSUB, NQS, IERR) 
         IF (IERR /= 0) GO TO 26 
!
!   Read the X, J, and (sign of) P quantum numbers
!
         READ (21, '(A)', IOSTAT=IOS) RECORD 
         IF (IOS /= 0) THEN 
            WRITE (6, *) 'LODCSL: Expecting intermediate' 
            WRITE (6, *) ' and final angular momentum' 
            WRITE (6, *) ' quantum number and final parity' 
            WRITE (6, *) ' specification;' 
            GO TO 26 
         ENDIF 
!
!   Allocate additional storage if necessary
!
         IF (NCF > NCFD) THEN 
            NEWSIZ = NCFD + NCFD/2 
            CALL RALLOC (IQA, NNNW, NEWSIZ, 'IQA', 'LODCSLM')
            CALL RALLOC (JQSA, NNNW, 3, NEWSIZ, 'JQSA', 'LODCSLM')
            CALL RALLOC (JCUPA, NNNW, NEWSIZ, 'JCUPA', 'LODCSLM')
            NCFD = NEWSIZ 
         ENDIF 
!
!   Zero out the arrays that store packed integers
!
         IQA(:NNNW,NCF) = 0 
         JQSA(:NNNW,1,NCF) = 0 
         JQSA(:NNNW,2,NCF) = 0 
         JQSA(:NNNW,3,NCF) = 0 
         JCUPA(:NNNW,NCF) = 0 
!
!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         DO I = 256, 1, -1 
            IF (RECORD(I:I) == ' ') CYCLE  
            LOC = I 
            EXIT  
         END DO 
         RECL = RECORD(LOC:LOC) 
         IF (RECL == '+') THEN 
            ISPARC = 1 
         ELSE IF (RECL == '-') THEN 
            ISPARC = -1 
         ELSE 
            WRITE (6, *) 'LODCSL: Incorrect parity' 
            WRITE (6, *) ' specification;' 
            GO TO 26 
         ENDIF 
         LOC = LOC - 1 
!
         CALL PARSJL (2, NCORE, RECORD, LOC, JX, NJX, IERR) 
         IF (IERR /= 0) GO TO 26 
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         DO I = 1, NCORE 
            CALL PACK (NKJ(I) + 1, I, IQA(1:NNNW,NCF)) 
            CALL PACK (0, I, JQSA(1:NNNW,1,NCF)) 
            CALL PACK (0, I, JQSA(1:NNNW,2,NCF)) 
            CALL PACK (1, I, JQSA(1:NNNW,3,NCF)) 
         END DO 
!
!   Check all subshell, intermediate and final angular momentum
!   quantum numbers; set the array elements in IQ, JQS for the peel
!   subshells; set the coupling array element in JCUP and the total
!   angular momentum array element in ITJPO
!
         IOC = 0 
         IPTY = 0 
         NQSN = 0 
         NJXN = 0 
         NPEELN = 0 
         NOPEN = 0 
         JLAST = 0 
         ILAST = 0 
         DO I = NCORP1, NW 
            IOCCI = IOCC(I) 
            NPEELN = NPEELN + IOCCI 
            NKJI = NKJ(I) 
            IFULLI = NKJI + 1 
            EMPTY = IOCCI == 0 
            IF (.NOT.EMPTY) IOC = IOC + 1 
            FULL = IOCCI == IFULLI 
            IF (EMPTY .OR. FULL) THEN 
               NU = 0 
               JSUB = 0 
            ELSE 
               IPTY = IPTY + NKL(I)*IOCCI 
               IF (NKJI /= 7) THEN 
                  NQSN = NQSN + 1 
                  IF (NQSN > NQS) THEN 
                     WRITE (6, *) 'LODCSL: Too few subshell' 
                     WRITE (6, *) ' quantum numbers specified;' 
                     GO TO 26 
                  ENDIF 
                  NU = 0 
                  JSUB = IQSUB(NQSN) 
               ELSE 
                  IF (IOCCI /= 4) THEN 
                     NQSN = NQSN + 1 
                     IF (NQSN > NQS) THEN 
                        WRITE (6, *) 'LODCSL: Too few subshell' 
                        WRITE (6, *) ' quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     NU = 0 
                     JSUB = IQSUB(NQSN) 
                  ELSE 
                     NQSN = NQSN + 1 
                     IF (NQSN > NQS) THEN 
                        WRITE (6, *) 'LODCSL: Too few subshell' 
                        WRITE (6, *) ' quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     JSUB = IQSUB(NQSN) 
                     IF (JSUB==4 .OR. JSUB==8) THEN 
                        NU = JSUB/2 
                        NQSN = NQSN + 1 
                        IF (NQSN > NQS) THEN 
                           WRITE (6, *) 'LODCSL: Too few subshell' 
                           WRITE (6, *) ' quantum numbers specified;' 
                           GO TO 26 
                        ENDIF 
                        JSUB = IQSUB(NQSN) 
                     ELSE 
                        NU = 0 
                     ENDIF 
                  ENDIF 
               ENDIF 
               IQT = MIN(IOCCI,IFULLI - IOCCI) 
               LOC = (IFULLI - 2)/2 
               LOC = (LOC*(LOC + 1))/2 + IQT 
               NBEG = JTAB(LOC+1) + 1 
               NEND = JTAB(LOC+2) 
               DO J = NBEG, NEND, 3 
                  IF (NTAB(J+2) /= JSUB + 1) CYCLE  
                  IF (NU == 0) THEN 
                     NU = NTAB(J) 
                     GO TO 9 
                  ELSE 
                     IF (NTAB(J) == NU) GO TO 9 
                  ENDIF 
               END DO 
               CALL CONVRT (NP(I), RECORD, LENTH) 
               WRITE (6, *) 'LODCSL: Subshell quantum numbers' 
               WRITE (6, *) ' specified incorrectly for' 
               WRITE (6, *) ' '//RECORD(1:LENTH)//NH(I)//' subshell.' 
               GO TO 26 
            ENDIF 
    9       CONTINUE 
            IF (.NOT.EMPTY .AND. .NOT.FULL) THEN 
               NOPEN = NOPEN + 1 
               IF (NOPEN > 1) THEN 
                  IF (JSUB == 0) THEN 
                     JXN = JLAST 
                  ELSE 
                     ILAST = IOC 
                     NJXN = NJXN + 1 
                     IF (NJXN > NJX) THEN 
                        WRITE (6, *) 'LODCSL: Too few intermediate' 
                        WRITE (6, *) ' and final angular momentum' 
                        WRITE (6, *) ' quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     JXN = JX(NJXN) 
                     DO J = ABS(JLAST - JSUB), JLAST + JSUB, 2 
                        IF (JXN == J) GO TO 11 
                     END DO 
                     CALL CONVRT (NP(I), RECORD, LENTH) 
                     WRITE (6, *) 'LODCSL: coupling of '//RECORD(1:LENTH)//NH(I&
                        ) 
                     WRITE (6, *) ' subshell to previous subshells' 
                     WRITE (6, *) ' is incorrect.' 
                     GO TO 26 
                  ENDIF 
   11             CONTINUE 
                  CALL PACK (JXN + 1, NOPEN - 1, JCUPA(1:NNNW,NCF)) 
                  JLAST = JXN 
               ELSE 
                  JLAST = JSUB 
               ENDIF 
            ENDIF 
            CALL PACK (IOCCI, I, IQA(1:NNNW,NCF)) 
            CALL PACK (NU, I, JQSA(1:NNNW,1,NCF)) 
            CALL PACK (0, I, JQSA(1:NNNW,2,NCF)) 
            CALL PACK (JSUB + 1, I, JQSA(1:NNNW,3,NCF)) 
         END DO 
!
         DO I = MAX(1,NOPEN), NW 
            CALL PACK (0, I, JCUPA(1:NNNW,NCF)) 
         END DO 
!
         IF (NQSN /= NQS) THEN 
            WRITE (6, *) 'LODCSL: Too many subshell' 
            WRITE (6, *) ' quantum numbers specified;' 
            GO TO 26 
         ENDIF 
!
         IF (ILAST /= IOC) NJXN = NJXN + 1 
         IF (NJXN /= NJX) THEN 
            WRITE (6, *) 'LODCSL: Too many intermediate' 
            WRITE (6, *) ' and final angular momentum' 
            WRITE (6, *) ' quantum numbers specified;' 
            GO TO 26 
         ENDIF 
!
         IF (JX(NJXN) /= JLAST) THEN 
            WRITE (6, *) 'LODCSL: Final angular momentum' 
            WRITE (6, *) ' incorrectly specified;' 
            GO TO 26 
         ENDIF 
!
         IPTY = (-1)**IPTY 
         IF (IPTY /= ISPARC) THEN 
            WRITE (6, *) 'LODCSL: Parity specified incorrectly;' 
            GO TO 26 
         ENDIF 
!
         JPI = (JLAST + 1)*IPTY 
         CALL PACK (JPI, NNNW, JCUPA(1:NNNW,NCF)) 
!
         IF (NCF > 1) THEN 
            IF (NPEELN /= NPEEL) THEN 
               WRITE (6, *) 'LODCSL: Inconsistency in the number' 
               WRITE (6, *) ' of electrons.' 
               GO TO 26 
            ENDIF 
         ELSE 
            NPEEL = NPEELN 
         ENDIF 
!
!   Check if this CSF was already in the list; stop with a
!   message if this is the case
!
!         print *, 'Check duplicated CSFs'
         IF (NCF > 1) THEN 
            DO J = 1, NCF - 1 
!                  print *,'j= ',j,ncf
               DO I = NCORP1, NW 
!                  print *,i
!                  print *, IQ(I,J), JQS(1,I,J), JQS(2,I,J),JQS(3,I,J)
!          print *, IQ(I,ncf), JQS(1,I,ncf), JQS(2,I,ncf),JQS(3,I,ncf)
                  IF (IQ(I,J) /= IQ(I,NCF)) GO TO 17 
                  IF (JQS(1,I,J) /= JQS(1,I,NCF)) GO TO 17 
                  IF (JQS(2,I,J) /= JQS(2,I,NCF)) GO TO 17 
                  IF (JQS(3,I,J) /= JQS(3,I,NCF)) GO TO 17 
               END DO 
               DO I = 1, NOPEN - 1 
                  WRITE (6, *) I 
                  WRITE (6, *) JCUP(I,J), JCUP(I,NCF) 
                  IF (JCUP(I,J) /= JCUP(I,NCF)) GO TO 17 
               END DO 
            END DO 
            WRITE (6, *) 'LODCSL: Repeated CSF;' 
            GO TO 26 
         ENDIF 
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    CONTINUE 
         NREC = NREC + 3 
         GO TO 3 
!
      ELSE 
!
!   There is always at least one CSF
!
         IF (NCF == 1) THEN 
            DO I = 1, NCORE 
               CALL PACK (NKJ(I) + 1, I, IQA(1:NNNW,1)) 
               CALL PACK (0, I, JQSA(1:NNNW,1,1)) 
               CALL PACK (0, I, JQSA(1:NNNW,2,1)) 
               CALL PACK (1, I, JQSA(1:NNNW,3,1)) 
            END DO 
            CALL PACK (0, 1, JCUPA(1:NNNW,1)) 
            CALL PACK (1, NNNW, JCUPA(1:NNNW,1)) 
         ELSE 
            NCF = NCF - 1 
         ENDIF 
!
      ENDIF 
!
!   Check if any subshell is empty; eliminate it from the
!   list if this is the case; issue a message
!
      I = NCORP1 
   19 CONTINUE 
      IF (I <= NW) THEN 
         DO J = 1, NCF 
            IF (IQ(I,J) /= 0) GO TO 23 
         END DO 
         CALL CONVRT (NP(I), RECORD, LENTH) 
         WRITE (6, *) 'Subshell '//RECORD(1:LENTH)//NH(I)//' is empty' 
         WRITE (6, *) ' in all CSFs; eliminating this' 
         WRITE (6, *) ' subshell from the list;' 
         NW = NW - 1 
         DO II = I, NW 
            NP(II) = NP(II+1) 
            NAK(II) = NAK(II+1) 
            NKL(II) = NKL(II+1) 
            NKJ(II) = NKJ(II+1) 
            NH(II) = NH(II+1) 
            DO J = 1, NCF 
               ITEMP = IQ(II + 1,J) 
               CALL PACK (ITEMP, II, IQA(1:NNNW,J)) 
               ITEMP = JQS(1,II + 1,J) 
               CALL PACK (ITEMP, II, JQSA(II:NNNW,1,J)) 
               ITEMP = JQS(2,II + 1,J) 
               CALL PACK (ITEMP, II, JQSA(II:NNNW,2,J)) 
               ITEMP = JQS(3,II + 1,J) 
               CALL PACK (ITEMP, II, JQSA(II:NNNW,3,J)) 
            END DO 
         END DO 
   23    CONTINUE 
         I = I + 1 
         GO TO 19 
      ENDIF 
!
!   Store the number of electrons in the COMMON variable
!
      NCOREL = 0 
      NCOREL = SUM(NKJ(:NCORE)+1) 
      NELEC = NCOREL + NPEEL 
!
!   All done; report
!
      CALL CONVRT (NCF, RECORD, LENTH) 
      WRITE (6, *) ' there are '//RECORD(1:LENTH)//' relativistic CSFs;' 
      WRITE (6, *) ' ... load complete;' 
!
!   Debug printout
!
      IF (LDBPA(1)) THEN 
         WRITE (99, *) 'From LODCSL:' 
         DO I = 1, NCF 
            WRITE (99, *) 'CSF ', I 
            WRITE (99, *) 'ITJPO: ', ITJPO(I) 
            WRITE (99, *) 'ISPAR: ', ISPAR(I) 
            WRITE (99, *) 'IQ: ', (IQ(J,I),J=1,NW) 
            WRITE (99, *) 'JQS(1): ', (JQS(1,J,I),J=1,NW) 
            WRITE (99, *) 'JQS(2): ', (JQS(2,J,I),J=1,NW) 
            WRITE (99, *) 'JQS(3): ', (JQS(3,J,I),J=1,NW) 
            WRITE (99, *) 'JCUP: ', (JCUP(J,I),J=1,NW - 1) 
         END DO 
      ENDIF 
!
      RETURN  
!
   26 CONTINUE 
      CALL CONVRT (NCF, RECORD, LENTH) 
      WRITE (6, *) ' CSF sequence number: '//RECORD(1:LENTH)//':' 
      REWIND (21) 
      DO I = 1, NREC 
         READ (21, *) 
      END DO 
      DO I = 1, 3 
         READ (21, '(A)', ERR=29, END=29) RECORD 
         LENTH = LEN_TRIM(RECORD) 
         WRITE (6, *) RECORD(1:LENTH) 
      END DO 
   29 CONTINUE 
      CLOSE(21) 
      STOP  
!
      END SUBROUTINE LODCSLM 
