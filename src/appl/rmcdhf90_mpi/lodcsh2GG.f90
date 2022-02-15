!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSH2GG(NFILE, NCORE, JB)
!
! IMPORTANT:
! ==========
!   If jb == LOADALL, then it loads ALL blocks. In this case, the
!   NCFblock in the common should be the total ncf and allocations
!   should be changed to this ncf. Both should be done before calling
!   this routine - like the genuine block case.
!
!   Loads the data from the  .csl  file for the current block          *
!   All parameters are inputs:
!
!     nfile:    The unit number, usually 21
!     jb:       The block number (see above)
!     ncore:    The number of the core (sub-)shells.
!   Since NCFblock is known in the block version, allocation is done
!   once outside this routine. NCFblock is checked against the one (NCF)
!   obtained here.
!
!   Call(s) to: [LIB92]: CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,
!                        PACK, PARSJL, PRSRCN, PRSRSL
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by C. F. Fischer for block calculation      22 May 1997   *
!   Updated by Xinghong He                               08 Jul 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:13:05   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE DEBUG_C
      USE DEF_C
      USE ORB_C, ncfblock => ncf
      USE SYMA_C,          ONLY: JPGG
      USE TERMS_C,         only: jtab, ntab
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE prsrcn_I
      USE parsjl_I
      USE pack_I
      USE convrt_I
      USE iq_I
      USE jqs_I
      USE jcup_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NFILE
      INTEGER  :: NCORE
      INTEGER, INTENT(IN) :: JB
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LOADALL = -119
      CHARACTER*7, PARAMETER :: MYNAME = 'LODCSH2'
      INTEGER, PARAMETER :: NW2 = 2*NNNW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NNNW) :: IOCC
      INTEGER, DIMENSION(NW2) :: IQSUB
      INTEGER, DIMENSION(NNNW) :: JX
      INTEGER :: NCORP1, NREC, NCF, NPEEL, I, J
      INTEGER :: IOS, IERR, LOC, NQS, ISPARC, NJX, IOC, IPTY
      INTEGER :: NQSN, NJXN, NPEELN, NOPEN, JLAST, ILAST, IOCCI
      INTEGER :: NKJI, IFULLI, NU, JSUB, IQT, NBEG, NEND
      INTEGER :: LENTH, JXN, JPI, NCOREL, IQGG, JBGG, NCFGG
      LOGICAL :: EMPTY, FULL
      CHARACTER :: STR*256, RECL
!-----------------------------------------------
!
      IF (JB /= LOADALL) THEN
         WRITE (6, *) 'Loading CSF File for block ', JB
      ELSE
         WRITE (6, *) 'Loading CSF File for ALL blocks '
      ENDIF

      NCORP1 = NCORE + 1
      NPEEL = NW - NCORE
!
! NPEEL is used as 1) number of peel orbitals (here) and
!                  2) number of peel electrons (later in this routine)
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
!
!   Zero out the arrays that store packed integers - only when ncfblock>0
!
      IQA(:NNNW,:NCFBLOCK) = 0
!GG      JQSA(:NNNW,1,:NCFBLOCK) = 0
!GG      JQSA(:NNNW,2,:NCFBLOCK) = 0
!GG      JQSA(:NNNW,3,:NCFBLOCK) = 0
!GG      JCUPA(:NNNW,:NCFBLOCK) = 0

      NCF = 0
!GGGG
      NCFGG = 0
      JBGG = 1
    3 NCF = NCF + 1
      NCFGG = NCFGG + 1
!GG    3 CONTINUE
!GG      NCF = NCF + 1
!GGGG
!
      READ (NFILE, '(A)', IOSTAT=IOS) STR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This IF...READ makes the routine load the entire file (all blocks)
! by ignoring the end-of-block mark

!GGGG
!GG      IF (IOS .EQ. 0 .AND. str(1:2) .EQ. ' *' .AND. jb .EQ. LOADALL)
!&
!GG        READ (nfile, '(A)', IOSTAT = IOS) str
      IF (IOS.EQ.0 .AND. str(1:2).EQ.' *' .AND. jb.EQ.LOADALL) THEN
         READ (nfile, '(A)', IOSTAT = IOS) str
         NCFGG = 1
         JBGG = JBGG + 1
      END IF
!GGGG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IOS==0 .AND. STR(1:2)/=' *') THEN
!
!   Read in the occupations (q) of the peel shells; stop with a
!   message if an error occurs
!
         CALL PRSRCN (STR, NCORE, IOCC, IERR)
         IF (IERR /= 0) GO TO 28
!
!   Read the J_sub and v quantum numbers
!
         READ (nfile,'(A)',IOSTAT = IOS) str
         IF (IOS /= 0) THEN
            WRITE (ISTDE, *) MYNAME//': Expecting subshell quantum', &
               ' number specification;'
            GO TO 27
         ENDIF
         LOC = LEN_TRIM(STR)
         CALL PARSJL (1, NCORE, STR, LOC, IQSUB, NQS, IERR)
         IF (IERR /= 0) GO TO 27
!
!   Read the X, J, and (sign of) P quantum numbers
!
         READ (nfile,'(A)',IOSTAT = IOS) str
         IF (IOS /= 0) THEN
            WRITE (ISTDE, *) MYNAME//': Expecting intermediate ', &
               'and final angular momentum'
            WRITE (ISTDE, *) 'quantum number and final parity ', &
               'specification;'
            GO TO 26
         ENDIF
!
!   Zero out the arrays that store packed integers
!
         IQA(:NNNW,NCF) = 0
!GG         JQSA(:NNNW,1,NCF) = 0
!GG         JQSA(:NNNW,2,NCF) = 0
!GG         JQSA(:NNNW,3,NCF) = 0
!GG         JCUPA(:NNNW,NCF) = 0

!   Determine the parity and all intermediate and the final
!   angular momentum quantum numbers
!
         LOC = LEN_TRIM(STR)
         RECL = STR(LOC:LOC)
         IF (RECL == '+') THEN
            ISPARC = 1
         ELSE IF (RECL == '-') THEN
            ISPARC = -1
         ELSE
            WRITE (ISTDE, *) MYNAME//': Incorrect parity ', &
                            'specification;'
            GO TO 26
         ENDIF
         LOC = LOC - 1
!
         CALL PARSJL (2, NCORE, STR, LOC, JX, NJX, IERR)
         IF (IERR /= 0) GO TO 26
!
!   Set the occupation and subshell quantum number array elements
!   in IQ, JQS for the core subshells
!
         DO I = 1, NCORE
            CALL PACK (NKJ(I) + 1, I, IQA(1:NNNW,NCF))
!GG            CALL PACK (0, I, JQSA(1:NNNW,1,NCF))
!GG            CALL PACK (0, I, JQSA(1:NNNW,2,NCF))
!GG            CALL PACK (1, I, JQSA(1:NNNW,3,NCF))
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
                     WRITE (ISTDE, *) MYNAME//': Too few subshell ', &
                        'quantum numbers specified;'
                     GO TO 26
                  ENDIF
                  NU = 0
                  JSUB = IQSUB(NQSN)
               ELSE
                  IF (IOCCI /= 4) THEN
                     NQSN = NQSN + 1
                     IF (NQSN > NQS) THEN
                        WRITE (ISTDE, *) MYNAME//': Too few subshell ', &
                           'quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     NU = 0
                     JSUB = IQSUB(NQSN)
                  ELSE
                     NQSN = NQSN + 1
                     IF (NQSN > NQS) THEN
                        WRITE (ISTDE, *) MYNAME//': Too few subshell ', &
                           'quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     JSUB = IQSUB(NQSN)
                     IF (JSUB==4 .OR. JSUB==8) THEN
                        NU = JSUB/2
                        NQSN = NQSN + 1
                        IF (NQSN > NQS) THEN
                           WRITE (ISTDE, *) MYNAME//': Too few subshell', &
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
               CALL CONVRT (NP(I), STR, LENTH)
               WRITE (ISTDE, *) MYNAME//': Subshell quantum numbers ', &
                  'specified incorrectly for '//STR(1:LENTH)//NH(I)//&
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
                        WRITE (ISTDE, *) MYNAME//': Too few intermediate', &
                           ' and final angular momentum', &
                           ' quantum numbers specified;'
                        GO TO 26
                     ENDIF
                     JXN = JX(NJXN)
                     DO J = ABS(JLAST - JSUB), JLAST + JSUB, 2
                        IF (JXN == J) GO TO 11
                     END DO
                     CALL CONVRT (NP(I), STR, LENTH)
                     WRITE (ISTDE, *) MYNAME//': coupling of '//STR(1:LENTH)//&
                        NH(I), ' subshell to previous subshells is incorrect.'
                     GO TO 26
                  ENDIF
   11             CONTINUE
!GG                  CALL PACK (JXN + 1, NOPEN - 1, JCUPA(1:NNNW,NCF))
                  JLAST = JXN
               ELSE
                  JLAST = JSUB
               ENDIF
            ENDIF
            CALL PACK(IOCCI, I, IQA(1:NNNW,NCF))
!GG            CALL PACK (NU, I, JQSA(1:NNNW,1,NCF))
!GG            CALL PACK (0, I, JQSA(1:NNNW,2,NCF))
!GG            CALL PACK (JSUB + 1, I, JQSA(1:NNNW,3,NCF))
         END DO
!
!GG         DO I = MAX(1,NOPEN), NW
!GG            CALL PACK (0, I, JCUPA(1:NNNW,NCF))
!GG         END DO
!
         IF (NQSN /= NQS) THEN
            WRITE (ISTDE, *) MYNAME//': Too many subshell', &
               ' quantum numbers specified;'
            GO TO 26
         ENDIF
!
         IF (ILAST /= IOC) NJXN = NJXN + 1
         IF (NJXN /= NJX) THEN
            WRITE (ISTDE, *) MYNAME//': Too many intermediate', &
               ' and final angular momentum', ' quantum numbers specified;'
            GO TO 26
         ENDIF
!
         IF (JX(NJXN) /= JLAST) THEN
            WRITE (ISTDE, *) MYNAME//': Final angular momentum', &
               ' incorrectly specified;'
            GO TO 26
         ENDIF
!
         IPTY = (-1)**IPTY
         IF (IPTY /= ISPARC) THEN
            WRITE (ISTDE, *) MYNAME//': Parity specified incorrectly;'
            GO TO 26
         ENDIF
!
         JPI = (JLAST + 1)*IPTY
!GGGG
         IF(NCFGG .EQ. 1) THEN
            JPGG(JBGG) = JPI
         END IF
!GG         CALL PACK (JPI, NNNW, JCUPA(1:NNNW,NCF))
!GGGG
!
         IF (NCF > 1) THEN
            IF (NPEELN /= NPEEL) THEN
               WRITE (ISTDE, *) MYNAME//': Inconsistency in the number', &
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
!GG         IF (NCF > 1) THEN
!GG            DO J = 1, NCF - 1
!GG               DO I = NCORP1, NW
!GG                  IF (IQ(I,J) /= IQ(I,NCF)) GO TO 17
!GG                  IF (JQS(1,I,J) /= JQS(1,I,NCF)) GO TO 17
!GG                  IF (JQS(2,I,J) /= JQS(2,I,NCF)) GO TO 17
!GG                  IF (JQS(3,I,J) /= JQS(3,I,NCF)) GO TO 17
!GG               END DO
!GG               DO I = 1, NOPEN - 1
!GG                  IF (JCUP(I,J) /= JCUP(I,NCF)) GO TO 17
!GG               END DO
!GG            END DO
!GG            WRITE (ISTDE, *) MYNAME//': Repeated CSF;'
!GG            GO TO 26
!GG         ENDIF
!
!   Successfully read a CSF; update NREC and read another CSF
!
   17    CONTINUE
         NREC = NREC + 3

         GO TO 3
!
      ELSE                                       ! the record just read is either ' *' or EOF, marking
            ! the end of a block or end of the file
!
!   There is always at least one CSF
!
         IF (NCF == 1) THEN
            DO I = 1, NCORE
               CALL PACK (NKJ(I) + 1, I, IQA(1:NNNW,1))
!GG               CALL PACK (0, I, JQSA(1:NNNW,1,1))
!GG               CALL PACK (0, I, JQSA(1:NNNW,2,1))
!GG               CALL PACK (1, I, JQSA(1:NNNW,3,1))
            END DO
!GG            CALL PACK (0, 1, JCUPA(1:NNNW,1))
!GG            CALL PACK (1, NNNW, JCUPA(1:NNNW,1))
         ELSE
            NCF = NCF - 1
         ENDIF
!
      ENDIF

      IF (NCF /= NCFBLOCK) THEN
         WRITE (ISTDE, *) MYNAME//': ncf=', NCF, 'ncfblock=', NCFBLOCK
         ERROR STOP
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
         CALL CONVRT (NP(I), STR, LENTH)
         WRITE (6, *) 'Subshell '//STR(1:LENTH)//NH(I)//' is empty', &
            ' in all CSFs'
   23    CONTINUE
         I = I + 1
         GO TO 19
      ENDIF
!
!   Store the number of electrons in the COMMON variable
!   This will act as a check now - it's been determined in lodcsh
!
      NCOREL = 0
      NCOREL = SUM(NKJ(:NCORE)+1)
!      NELEC = NCOREL+NPEEL
      IF (NCOREL + NPEEL /= NELEC) THEN
         WRITE (ISTDE, *) MYNAME//': nelec not equal to that in lodcsh'
         ERROR STOP
      ENDIF
      WRITE (6,*)'There are ',NCF,' relativistic CSFs... load complete;'
      RETURN
!
   26 CONTINUE
      BACKSPACE (NFILE)
   27 CONTINUE
      BACKSPACE (NFILE)
   28 CONTINUE
      BACKSPACE (NFILE)
      WRITE (ISTDE, *) ' CSF sequence number: ', NCF
      DO I = 1, 3
         READ  (nfile,'(A)',ERR = 29,END = 29) str
         WRITE (ISTDE, *) STR(1:LEN_TRIM(STR))
      END DO
   29 continue
      CLOSE(NFILE)

      ERROR STOP
      END SUBROUTINE LODCSH2GG
