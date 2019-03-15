!***********************************************************************
!                                                                      *
      SUBROUTINE SORTMEM(NFILE,IGGMAX_K,NCOEFF,NTGRAL,LPRINT,NB,FHEAD)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE memory_man
      USE ORB_C,           ONLY: NP, NCF, NH
      USE IOUNIT_C
      USE sacoef_C, COEFF=>TCOEFF,LABEL=>ILABEL,ICLMN=>IICLMN,INDEX=>IINDEX
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
      INTEGER, INTENT(IN)   :: NB,IGGMAX_K
      LOGICAL, INTENT(IN)   :: LPRINT
      CHARACTER, INTENT(IN) :: FHEAD*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, LABEL, NSWAP
      INTEGER, DIMENSION(:), pointer :: NSWAP
!      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
      INTEGER :: LCNUM, LCK, IERR, MB, I, L, IR, ICL, IND, LAB, NSW, J,&
         LAST, IBEG, IEND, NCONTR, IA, IB, K, ID, IC, IGG
      REAL(DOUBLE) :: COF
      CHARACTER :: CNUM*20, SRTLAB*8, MCPLAB*3, CK*2
      CHARACTER (LEN = LEN (fhead) + 3):: fullname
!-----------------------------------------------
!
      IF(NFILE .EQ. 31) THEN
         IGG = 1
      ELSE
         IGG = NFILE - 30
      END IF
!      IGG = IGGMAX_K*(IGG-1)
!   Message
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

      FULLNAME = FHEAD//'.'//CK(1:2)

      OPEN(29,FILE=FULLNAME,STATUS='OLD',FORM='UNFORMATTED', &
         IOSTAT=IERR, POSITION='APPEND')
      IF (IERR /= 0) THEN
         WRITE (ISTDE, *) ' Error when opening the file ', FULLNAME
         STOP
      ENDIF
!
!   Sort the list
!
      IF (NCOEFF > 0) THEN
!
!   Allocate storage for all required arrays
!
         IF (NFILE == 33) CALL ALLOC (NSWAP, NCOEFF,'NSWAP', 'SORTMEM')
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
            COF = COEFF(L,IGG)
            ICL = ICLMN(L,IGG)
            IND = INDEX(L,IGG)
            LAB = LABEL(L,IGG)
            NSW = NSWAP(L)
         ELSE
            COF = COEFF(IR,IGG)
            ICL = ICLMN(IR,IGG)
            IND = INDEX(IR,IGG)
            LAB = LABEL(IR,IGG)
            NSW = NSWAP(IR)
            COEFF(IR,IGG) = COEFF(1,IGG)
            ICLMN(IR,IGG) = ICLMN(1,IGG)
            INDEX(IR,IGG) = INDEX(1,IGG)
            LABEL(IR,IGG) = LABEL(1,IGG)
            NSWAP(IR) = NSWAP(1)
            IR = IR - 1
            IF (IR == 1) THEN
               COEFF(1,IGG) = COF
               ICLMN(1,IGG) = ICL
               INDEX(1,IGG) = IND
               LABEL(1,IGG) = LAB
               NSWAP(1) = NSW
               GO TO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
  345    CONTINUE
         IF (J <= IR) THEN
            IF (J < IR) THEN
               IF (LABEL(J,IGG) < LABEL(J+1,IGG)) J = J + 1
            ENDIF
            IF (LAB < LABEL(J,IGG)) THEN
               COEFF(I,IGG) = COEFF(J,IGG)
               ICLMN(I,IGG) = ICLMN(J,IGG)
               INDEX(I,IGG) = INDEX(J,IGG)
               LABEL(I,IGG) = LABEL(J,IGG)
               NSWAP(I) = NSWAP(J)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GO TO 345
         ENDIF
         COEFF(I,IGG) = COF
         ICLMN(I,IGG) = ICL
         INDEX(I,IGG) = IND
         LABEL(I,IGG) = LAB
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
            COF = COEFF(L,IGG)
            ICL = ICLMN(L,IGG)
            IND = INDEX(L,IGG)
            LAB = LABEL(L,IGG)
         ELSE
            COF = COEFF(IR,IGG)
            ICL = ICLMN(IR,IGG)
            IND = INDEX(IR,IGG)
            LAB = LABEL(IR,IGG)
            COEFF(IR,IGG) = COEFF(1,IGG)
            ICLMN(IR,IGG) = ICLMN(1,IGG)
            INDEX(IR,IGG) = INDEX(1,IGG)
            LABEL(IR,IGG) = LABEL(1,IGG)
            IR = IR - 1
            IF (IR == 1) THEN
               COEFF(1,IGG) = COF
               ICLMN(1,IGG) = ICL
               INDEX(1,IGG) = IND
               LABEL(1,IGG) = LAB
               GO TO 456
            ENDIF
         ENDIF
         I = L
         J = L + L
   93    CONTINUE
         IF (J <= IR) THEN
            IF (J < IR) THEN
               IF (LABEL(J,IGG) < LABEL(J+1,IGG)) J = J + 1
            ENDIF
            IF (LAB < LABEL(J,IGG)) THEN
               COEFF(I,IGG) = COEFF(J,IGG)
               ICLMN(I,IGG) = ICLMN(J,IGG)
               INDEX(I,IGG) = INDEX(J,IGG)
               LABEL(I,IGG) = LABEL(J,IGG)
               I = J
               J = J + J
            ELSE
               J = IR + 1
            ENDIF
            GO TO 93
         ENDIF
         COEFF(I,IGG) = COF
         ICLMN(I,IGG) = ICL
         INDEX(I,IGG) = IND
         LABEL(I,IGG) = LAB
         GO TO 92

      ENDIF
!
!   Sorting complete; rewrite the file header
!
!
  456 CONTINUE
      WRITE (29) 'MCP', NB, NCF, NCOEFF
!GG      WRITE (9999,*) 'MCP', NB, NCF, NCOEFF, NFILE
!
!   Write the sorted list to mcp.xx
!
      IF (NCOEFF > 0) THEN
!
         LAST = LABEL(1,IGG)
         IBEG = 1
         IEND = 1
         NTGRAL = 1
!
         DO I = 2, NCOEFF
            IF (LABEL(I,IGG) == LAST) THEN
               IEND = IEND + 1
            ELSE
               WRITE (29) LAST, IEND - IBEG + 1
!GG               WRITE (9999,*) LAST, IEND - IBEG + 1
               WRITE (29) (ICLMN(J,IGG),INDEX(J,IGG),COEFF(J,IGG),     &
                                                          J=IBEG,IEND)
!cjb-GG              DO J = IBEG,IEND
!GG               WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!cjb-GG              END DO

!              IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,IEND)
               NTGRAL = NTGRAL + 1
               LAST = LABEL(I,IGG)
               IBEG = IEND + 1
               IEND = IBEG
            ENDIF
         END DO
!
         IF (IBEG <= NCOEFF) THEN
            WRITE (29) LAST, NCOEFF - IBEG + 1
!GG            WRITE (9999,*) LAST, NCOEFF - IBEG + 1
            WRITE (29) (ICLMN(J,IGG),INDEX(J,IGG),COEFF(J,IGG),       &
                                                       J=IBEG,NCOEFF)
!GG          DO J=IBEG,NCOEFF
!GG            WRITE (9999,'(2I12,E25.15)') ICLMN(J),INDEX(J),COEFF(J)
!GG          END DO

!           IF (NFILE.EQ.33) WRITE (20) (NSWAP(J),J = IBEG,NCOEFF)
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
      !PRINT *, ' ... sort complete; ', ntgral, ' integrals;'
!
!   Debug printout
!
      IF (LPRINT) THEN
         WRITE (99, 300)
         WRITE (6, 300)
         IF (NCOEFF > 0) THEN
!
            LAST = LABEL(1,IGG)
            IBEG = 1
            IEND = 1
!
            DO I = 2, NCOEFF
               IF (LABEL(I,IGG) == LAST) IEND = IEND + 1
  567          CONTINUE
               IF (LABEL(I,IGG)==LAST .AND. I/=NCOEFF) CYCLE
               LAB = LAST
               NCONTR = IEND - IBEG + 1
               IF (NFILE == 31) THEN
                  IA = MOD(LAB,KEY)
                  IB = LAB/KEY
                  WRITE (99, 301) NP(IA), NH(IA), NP(IB), NH(IB)
                  DO J = IBEG, IEND
                     WRITE (99, 302) ICLMN(J,IGG),INDEX(J,IGG),       &
                                                        COEFF(J,IGG)
                  END DO
               ELSE
                  K = NFILE - 32
                  ID = MOD(LAB,KEY)
                  LAB = LAB/KEY
                  IB = MOD(LAB,KEY)
                  LAB = LAB/KEY
                  IC = MOD(LAB,KEY)
                  IA = LAB/KEY
                  WRITE (99, 304) K, NP(IA), NH(IA), NP(IB), NH(IB), NP(IC), NH&
                     (IC), NP(ID), NH(ID)
                  DO J = IBEG, IEND
                     WRITE (99, 305) K,ICLMN(J,IGG),INDEX(J,IGG),     &
                                                          COEFF(J,IGG)
                  END DO
               ENDIF
               LAST = LABEL(I,IGG)
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
         IF (NFILE == 33) CALL DALLOC (NSWAP, 'NSWAP', 'SORTMEM')
      ENDIF

  300 FORMAT(/,'From SORT:')
  301 FORMAT(' I(',1I2,1A2,',',1I2,1A2,'):')
  302 FORMAT('  T_[',1I2,',',1I4,'] = ',1P,D19.12)
  303 FORMAT('  Number of integrals is ',1I4)
  304 FORMAT(' R^[(',1I2,')] (',1I2,1A2,',',1I2,1A2,';',1I2,1A2,',',1I2,1A2,&
         '):')
  305 FORMAT('  V^[(',1I2,')]_[',1I8,',',1I8,'] = ',1P,D19.12)

      RETURN
      END SUBROUTINE SORTMEM
