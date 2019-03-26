!***********************************************************************
!                                                                      *
      SUBROUTINE SORT(NFILE, NCOEFF, NTGRAL, LPRINT, NB, FHEAD)
!                                                                      *
!   This routine sorts lists                                           *
!                                                                      *
!                     (ICLMN,INDEX,LABEL,COEFF)                        *
!   into lists                                                         *
!                     (LABEL,ICLMN,INDEX,COEFF)                        *
!                                                                      *
!   using Heapsort. File NFILE is closed by this routine.  NCOEFF is   *
!   the number of triads (INDEX, ...) and is an input. NTGRAL is the   *
!   number of different values of LABEL, and is an output. If LPRINT   *
!   is  .TRUE. , the contents of the sorted file are interpreted and   *
!   printed to unit 99.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:24:58   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE memory_man
      USE ORB_C,           ONLY: NP, NCF, NH
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NFILE
      INTEGER  :: NCOEFF
      INTEGER, INTENT(OUT)  :: NTGRAL
      INTEGER, INTENT(IN)   :: NB
      LOGICAL, INTENT(IN)   :: LPRINT
      CHARACTER, INTENT(IN) :: FHEAD*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, LABEL, NSWAP
      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
      INTEGER :: LCNUM, LCK, IERR, MB, I, L, IR, ICL, IND, LAB, NSW, J,&
         LAST, IBEG, IEND, NCONTR, IA, IB, K, ID, IC
      REAL(DOUBLE) :: COF
      CHARACTER :: CNUM*20, SRTLAB*8, MCPLAB*3, CK*2
      CHARACTER (LEN = LEN (fhead) + 3):: fullname
!-----------------------------------------------
!
      CALL CONVRT (NCOEFF, CNUM, LCNUM)
      IF (NFILE > 31) THEN
         CALL CONVRT (NFILE - 32, CK, LCK)
         WRITE (6, *) 'Sorting '//CNUM(1:LCNUM)//' V(k='//CK(1:LCK)//&
            ') coefficients ...', NFILE
      ELSE
         WRITE (6, *) 'Sorting '//CNUM(1:LCNUM)//' T coefficients ...', NFILE
      ENDIF
!
      CALL CONVRT (NFILE, CK, LCK)
      IF (LCK > 2) THEN
         WRITE (ISTDE, *) 'sort: nfile > 99; check fullname'
         STOP
      ENDIF
!
      FULLNAME = FHEAD//'.'//CK(1:2)
      OPEN(29,FILE=FULLNAME,STATUS='OLD',FORM='UNFORMATTED', &
         IOSTAT=IERR, POSITION='APPEND')
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) ' Error when opening the file ', FULLNAME
         STOP
      ENDIF
!
      READ (NFILE) MCPLAB, MB
      IF (NB /= MB) THEN
         WRITE (ISTDE, *) 'sort: nb = ', NB, '.NE. mb (=', MB, ')'
         STOP
      ENDIF
!
!   Sort the list
!
      IF (NCOEFF > 0) THEN
!
!   Allocate storage for all required arrays
!
         CALL ALLOC (COEFF, NCOEFF, 'COEFF', 'SORT')
         CALL ALLOC (ICLMN, NCOEFF, 'ICLMN', 'SORT')
         CALL ALLOC (INDEX, NCOEFF, 'INDEX', 'SORT')
         CALL ALLOC (LABEL, NCOEFF, 'LABEL', 'SORT')
         IF (NFILE == 33) CALL ALLOC (NSWAP, NCOEFF,'NSWAP', 'SORT'  )
!
!   Read arrays into memory from NFILE
!
         DO I = 1, NCOEFF
            READ (NFILE) ICLMN(I), INDEX(I), LABEL(I), COEFF(I)
         END DO
      ENDIF
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
      IF (NFILE==33 .AND. NCOEFF>1) THEN

         L = NCOEFF/2 + 1
         IR = NCOEFF
  234    CONTINUE
         IF (L > 1) THEN
            L = L - 1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
            NSW = NSWAP(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            NSW = NSWAP(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            NSWAP(IR) = NSWAP(1)
            IR = IR - 1
            IF (IR == 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               NSWAP(1) = NSW
               GO TO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
  345    CONTINUE
         IF (J <= IR) THEN
            IF (J < IR) THEN
               IF (LABEL(J) < LABEL(J+1)) J = J + 1
            ENDIF
            IF (LAB < LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               NSWAP(I) = NSWAP(J)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GO TO 345
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         NSWAP(I) = NSW
         GO TO 234

      ELSE IF (NFILE/=33 .AND. NCOEFF>1) THEN
!
!   Sort LABEL into ascending order using the heapsort algorithm;
!   move the associated members of COEFF and INDEX in the same
!   manner; the code below is adapted from Press et al.
!
         L = NCOEFF/2 + 1
         IR = NCOEFF
   92    CONTINUE
         IF (L > 1) THEN
            L = L - 1
            COF = COEFF(L)
            ICL = ICLMN(L)
            IND = INDEX(L)
            LAB = LABEL(L)
         ELSE
            COF = COEFF(IR)
            ICL = ICLMN(IR)
            IND = INDEX(IR)
            LAB = LABEL(IR)
            COEFF(IR) = COEFF(1)
            ICLMN(IR) = ICLMN(1)
            INDEX(IR) = INDEX(1)
            LABEL(IR) = LABEL(1)
            IR = IR - 1
            IF (IR == 1) THEN
               COEFF(1) = COF
               ICLMN(1) = ICL
               INDEX(1) = IND
               LABEL(1) = LAB
               GO TO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
   93    CONTINUE
         IF (J <= IR) THEN
            IF (J < IR) THEN
               IF (LABEL(J) < LABEL(J+1)) J = J + 1
            ENDIF
            IF (LAB < LABEL(J)) THEN
               COEFF(I) = COEFF(J)
               ICLMN(I) = ICLMN(J)
               INDEX(I) = INDEX(J)
               LABEL(I) = LABEL(J)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GO TO 93
         ENDIF
         COEFF(I) = COF
         ICLMN(I) = ICL
         INDEX(I) = IND
         LABEL(I) = LAB
         GO TO 92

      ENDIF
!
!   Sorting complete; rewrite the file header
!
  456 CONTINUE
      WRITE (29) 'MCP', NB, NCF, NCOEFF
!GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
!
!   Write the sorted list to mcp.xx
!
      IF (NCOEFF > 0) THEN
!
         LAST = LABEL(1)
         IBEG = 1
         IEND = 1
         NTGRAL = 1
!
         DO I = 2, NCOEFF
            IF (LABEL(I) == LAST) THEN
               IEND = IEND + 1
            ELSE
               WRITE (29) LAST, IEND - IBEG + 1
!GG               WRITE (9999,*) LAST, IEND - IBEG + 1
               WRITE (29) (ICLMN(J),INDEX(J),COEFF(J),J=IBEG,IEND)
!cjb-GG       DO J = IBEG,IEND
!GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!cjb-GG        END DO

               NTGRAL = NTGRAL + 1
               LAST = LABEL(I)
               IBEG = IEND + 1
               IEND = IBEG
            ENDIF
         END DO
!
         IF (IBEG <= NCOEFF) THEN
            WRITE (29) LAST, NCOEFF - IBEG + 1
!GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            WRITE (29) (ICLMN(J),INDEX(J),COEFF(J),J=IBEG,NCOEFF)
          DO J=IBEG,NCOEFF
!GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
          END DO

         ENDIF
!
      ELSE
!
         NTGRAL = 0
!
      ENDIF
!
!     write the terminator record for this block
!
      WRITE (29) 0, 0
!GG      WRITE (9999,*) 0, 0
      CLOSE(29)
!
!   Completion message
!
!      PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
!
!   Debug printout
!
      IF (LPRINT) THEN
         WRITE (99, 300)
         WRITE (6, 300)
         IF (NCOEFF > 0) THEN
!
            LAST = LABEL(1)
            IBEG = 1
            IEND = 1
!
            DO I = 2, NCOEFF
               IF (LABEL(I) == LAST) IEND = IEND + 1
  567          CONTINUE
               IF (LABEL(I)==LAST .AND. I/=NCOEFF) CYCLE
               LAB = LAST
               NCONTR = IEND - IBEG + 1
               IF (NFILE == 31) THEN
                  IA = MOD(LAB,KEY)
                  IB = LAB/KEY
                  WRITE (99, 301) NP(IA), NH(IA), NP(IB), NH(IB)
                  DO J = IBEG, IEND
                     WRITE (99, 302) ICLMN(J), INDEX(J), COEFF(J)
                  END DO
               ELSE
                  K = NFILE - 32
                  ID = MOD(LAB,KEY)
                  LAB = LAB/KEY
                  IB = MOD(LAB,KEY)
                  LAB = LAB/KEY
                  IC = MOD(LAB,KEY)
                  IA = LAB/KEY
                  WRITE (99, 304)K,NP(IA),NH(IA),NP(IB),NH(IB),NP(IC), &
                                 NH(IC),NP(ID),NH(ID)
                  DO J = IBEG, IEND
                     WRITE (99, 305)K,ICLMN(J),INDEX(J),COEFF(J)
                  END DO
               ENDIF
               LAST = LABEL(I)
               IBEG = IEND + 1
               IEND = IBEG
               IF (IEND == NCOEFF) GO TO 567
            END DO
         ENDIF
         WRITE (99, 303) NTGRAL
      ENDIF
!
!   Deallocate storage
!
      IF (NCOEFF > 0) THEN
         CALL DALLOC (COEFF,'COEFF', 'SORT')
         CALL DALLOC (ICLMN, 'ICLMN', 'SORT')
         CALL DALLOC (INDEX, 'INDEX', 'SORT')
         CALL DALLOC (LABEL, 'LABEL', 'SORT')
         IF (NFILE == 33) CALL DALLOC (NSWAP, 'NSWAP', 'SORT')
      ENDIF

  300 FORMAT(/,'From SORT:')
  301 FORMAT(' I(',1I2,1A2,',',1I2,1A2,'):')
  302 FORMAT('  T_[',1I2,',',1I4,'] = ',1P,D19.12)
  303 FORMAT('  Number of integrals is ',1I4)
  304 FORMAT(' R^[(',1I2,')] (',1I2,1A2,',',1I2,1A2,';',1I2,1A2,',',1I2,1A2,&
         '):')
  305 FORMAT('  V^[(',1I2,')]_[',1I8,',',1I8,'] = ',1P,D19.12)

      RETURN
      END SUBROUTINE SORT
