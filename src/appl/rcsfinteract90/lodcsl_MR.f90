!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSL_MR(NCORE,NPEEL,NCFD,NEXT_BLOCK) 
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
!                        PACK, PARSJL, PRSRCN, PRSRSL                  *
!                                                                      *
!   Written by  G. Gaigalas                      NIST, December 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE parameter_def,   ONLY:  NNNW
      USE DEBUG_C 
      USE DEF_C 
      USE ORB_C 
      USE STAT_C
      USE TERMS_C,          only: jtab, ntab 
      USE IOUNIT_C 
      USE BLK_C,            only: NBLOCK,NCFBLK
      USE memory_man
      USE rang_Int_C
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
      INTEGER,  INTENT(IN)    :: NCORE
      INTEGER,  INTENT(INOUT) :: NPEEL
      INTEGER,  INTENT(OUT)   :: NCFD
      LOGICAL,  INTENT(OUT)   :: NEXT_BLOCK 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NW2 = 2*NNNW 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NNNW) :: IOCC 
      INTEGER, DIMENSION(NW2)  :: IQSUB 
      INTEGER, DIMENSION(NNNW) :: JX 
      INTEGER :: I
      INTEGER :: NCORP1, J, NPJ, NAKJ, LENTH, NREC &
         , IOS, IERR, LOC, NQS, NEWSIZ, ISPARC, NJX, IOC, IPTY, NQSN   &
         , NJXN, NPEELN, NOPEN, JLAST, ILAST, IOCCI, NKJI, IFULLI, NU  &
         , JSUB, IQT, NBEG, NEND, JXN, JPI, II, ITEMP, NCOREL 
      LOGICAL :: EMPTY, FULL 
      CHARACTER          :: RECL 
      CHARACTER(LEN=256) :: RECORD 
!-----------------------------------------------
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      NEXT_BLOCK = .TRUE.
      NCFD = NUM_in_BLK(NBLOCK+1)
      allocate (Found(NCFD+1))
      allocate (C_shell(NCFD+1))
      allocate (C_quant(NCFD+1))
      allocate (C_coupl(NCFD+1))
      allocate (IQA(NNNW,NCFD+1))
      allocate (JQSA(NNNW,3,NCFD+1))
      allocate (JCUPA(NNNW,NCFD+1))
      NREC = 5 
      NCF = 0 
    3 CONTINUE 
      NCF = NCF + 1 
!
      READ (21, '(A)', IOSTAT=IOS) RECORD 
      IF (RECORD(1:2) == ' *') THEN 
         NBLOCK = NBLOCK + 1 
         NCFBLK(NBLOCK) = NCF - 1 
         NotFound = NCFBLK(NBLOCK)
         Found(1:NotFound) = 0
         NEXT_BLOCK = .TRUE.
         RETURN
      ENDIF 
      C_shell(NCF) = RECORD
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
         IF (IOS /= 0) THEN 
            WRITE (ISTDE, *) 'LODCSL: Expecting subshell quantum', &
               ' number specification;' 
            GO TO 26 
         ENDIF 
         READ (21, '(A)', IOSTAT=IOS) RECORD 
         C_quant(NCF) = RECORD
         LOC = LEN_TRIM(RECORD) 
         CALL PARSJL (1, NCORE, RECORD, LOC, IQSUB, NQS, IERR) 
         IF (IERR /= 0) GO TO 26 
!
!   Read the X, J, and (sign of) P quantum numbers
!
         READ (21, '(A)', IOSTAT=IOS) RECORD 
         IF (IOS /= 0) THEN 
            WRITE (ISTDE, *) 'LODCSL: Expecting intermediate ', &
               'and final angular momentum' 
            WRITE (ISTDE, *) 'quantum number and final parity ', &
               'specification;' 
            GO TO 26 
         ENDIF 
         C_coupl(NCF) = RECORD
         WRITE(22,'(A)') TRIM(C_shell(NCF))
         WRITE(22,'(A)') TRIM(C_quant(NCF))
         WRITE(22,'(A)') TRIM(C_coupl(NCF))
!
!   Zero out the arrays that store packed integers
!
         DO I = 1,NNNW
            IQA(I,NCF)    = 0
            JQSA(I,1,NCF) = 0
            JQSA(I,2,NCF) = 0
            JQSA(I,3,NCF) = 0
            JCUPA(I,NCF)  = 0
         END DO
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
            WRITE (ISTDE, *) 'LODCSL: Incorrect parity ', &
                             'specification;' 
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
         NCORP1 = NCORE + 1
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
                     WRITE (ISTDE, *) 'LODCSL: Too few subshell quantum', &
                        ' numbers specified;' 
                     GO TO 26 
                  ENDIF 
                  NU = 0 
                  JSUB = IQSUB(NQSN) 
               ELSE 
                  IF (IOCCI /= 4) THEN 
                     NQSN = NQSN + 1 
                     IF (NQSN > NQS) THEN 
                        WRITE (ISTDE, *) 'LODCSL: Too few subshell ', &
                           'quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     NU = 0 
                     JSUB = IQSUB(NQSN) 
                  ELSE 
                     NQSN = NQSN + 1 
                     IF (NQSN > NQS) THEN 
                        WRITE (ISTDE, *) 'LODCSL: Too few subshell ', &
                           'quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     JSUB = IQSUB(NQSN) 
                     IF (JSUB==4 .OR. JSUB==8) THEN 
                        NU = JSUB/2 
                        NQSN = NQSN + 1 
                        IF (NQSN > NQS) THEN 
                           WRITE (ISTDE, *) 'LODCSL: Too few subshell', &
                              ' quantum numbers specified;' 
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
               WRITE (ISTDE, *) 'LODCSL: Subshell quantum numbers ', &
                  'specified incorrectly for '//RECORD(1:LENTH)//NH(I)//&
                  ' subshell.' 
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
                        WRITE (ISTDE, *) 'LODCSL: Too few intermediate', &
                           ' and final angular momentum', &
                           ' quantum numbers specified;' 
                        GO TO 26 
                     ENDIF 
                     JXN = JX(NJXN) 
                     DO J = ABS(JLAST - JSUB), JLAST + JSUB, 2 
                        IF (JXN == J) GO TO 11 
                     END DO 
                     CALL CONVRT (NP(I), RECORD, LENTH) 
                     WRITE (ISTDE, *) &
                        'LODCSL: coupling of '//RECORD(1:LENTH)//NH(I),&
                        ' subshell to previous subshells is incorrect.' 
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
            WRITE (ISTDE, *) 'LODCSL: Too many subshell', &
               ' quantum numbers specified;' 
            GO TO 26 
         ENDIF 
!
         IF (ILAST /= IOC) NJXN = NJXN + 1 
         IF (NJXN /= NJX) THEN 
            WRITE (ISTDE, *) 'LODCSL: Too many intermediate', &
            ' and final angular momentum', ' quantum numbers specified;' 
            GO TO 26 
         ENDIF 
!
         IF (JX(NJXN) /= JLAST) THEN 
            WRITE (ISTDE, *) 'LODCSL: Final angular momentum', &
               ' incorrectly specified;' 
            GO TO 26 
         ENDIF 
!
         IPTY = (-1)**IPTY 
         IF (IPTY /= ISPARC) THEN 
            WRITE (ISTDE, *) 'LODCSL: Parity specified incorrectly;' 
            GO TO 26 
         ENDIF 
!
         JPI = (JLAST + 1)*IPTY 
         CALL PACK (JPI, NNNW, JCUPA(1:NNNW,NCF)) 
!
         IF (NCF > 1) THEN 
            IF (NPEELN /= NPEEL) THEN 
               WRITE (ISTDE, *) 'LODCSL: Inconsistency in the number', &
                  ' of electrons.' 
               GO TO 26 
            ENDIF 
         ELSE 
            NPEEL = NPEELN 
         ENDIF 
!
!   Check if this CSF was already in the list; stop with a
!   message if this is the case
!
         IF (NCF > 1) THEN 
            DO J = 1, NCF - 1 
               DO I = NCORP1, NW 
                  IF (IQ(I,J) /= IQ(I,NCF)) GO TO 17 
                  IF (JQS(1,I,J) /= JQS(1,I,NCF)) GO TO 17 
                  IF (JQS(2,I,J) /= JQS(2,I,NCF)) GO TO 17 
                  IF (JQS(3,I,J) /= JQS(3,I,NCF)) GO TO 17 
               END DO 
               DO I = 1, NOPEN - 1 
                  IF (JCUP(I,J) /= JCUP(I,NCF)) GO TO 17 
               END DO 
            END DO 
            WRITE (ISTDE, *) 'LODCSL: Repeated CSF;' 
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
      I = NCORP1 
!
!   Store the number of electrons in the COMMON variable
!
      NCOREL = 0 
      NCOREL = SUM(NKJ(:NCORE)+1) 
      NELEC = NCOREL + NPEEL 
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
      NBLOCK = NBLOCK + 1 
      NCFBLK(NBLOCK) = NCF 
      NotFound = NCFBLK(NBLOCK)
      Found(1:NotFound) = 0
      NEXT_BLOCK = .FALSE.
!
      RETURN  
!
   26 CONTINUE 
      CALL CONVRT (NCF, RECORD, LENTH) 
      WRITE (ISTDE, *) ' CSF sequence number: '//RECORD(1:LENTH)//':' 
      REWIND (21) 
      DO I = 1, NREC 
         READ (21, *) 
      END DO 
      DO I = 1, 3 
         READ (21,'(A)',ERR = 29,END = 29) RECORD
         LENTH = LEN_TRIM(RECORD) 
         WRITE (ISTDE, *) RECORD(1:LENTH) 
      END DO 
   29 CLOSE(21) 
      STOP  
!
      END SUBROUTINE LODCSL_MR
