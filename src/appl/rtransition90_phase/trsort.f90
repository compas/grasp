!***********************************************************************
!                                                                      *
      SUBROUTINE TRSORT(NAME, NFILE, NFILE2, LPRINT, JKP, IBLKI, IBLKF)
!                                                                      *
!   Routine to sort angular coefficients into list based on integral   *
!   labels rather than CSF.  A tree sort is used. To save space, the   *
!   CSF pair labels and coefficients  are not read into arrays until   *
!   the sorting has been done.                                         *
!                                                                      *
!   Call(s) to: [OSCL92]: ALCLLA, ALCNMA.                              *
!                                                                      *
!   Original author(s) unknown. Modifications for dynamic storage by   *
!   Farid A. Parpia.                                                   *
!                                                                      *
!                                           Last update: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:53:13   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE memory_man
      USE default_C
      USE Orb_C
      USE osc_C, ONLY: nkp, kp
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alclla_I
      USE alcnma_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NFILE
      INTEGER, INTENT(IN) :: NFILE2
      INTEGER  :: JKP
      INTEGER, INTENT(IN) :: IBLKI
      INTEGER, INTENT(IN) :: IBLKF
      LOGICAL, INTENT(IN) :: LPRINT
      CHARACTER, INTENT(IN) :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NCA = 65536
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(:), pointer :: ILAB, IRIGHT, ILEFT, IBEG, IPTR, &
                    ILAST, IBLINT, IPTCSF, LBLINT, JLABL, ISLDR, ISLDR1
      REAL(DOUBLE), DIMENSION(:), pointer :: XSLDR, XL
      INTEGER :: NMCP, NINT, NLABEL, LLDIM, NMDIM, IR, IS, NI, I, M, J, &
                 ICOUNT, JLAB, K, L, IMCP, J1, J2, MLR, MUP, INTS, IA, IB
      REAL(DOUBLE) :: X
      LOGICAL :: FIRST
      CHARACTER(LEN=2), DIMENSION(-9:9) :: S
!-----------------------------------------------
!
!
!     POINTER (PIBEG,IBEG(1)),(PILAB,ILAB(1)),(PILAST,ILAST(1)),
!    :        (PILEFT,ILEFT(1)),(PIPTCS,IPTCSF(1)),(PIRIGH,IRIGHT(1)),
!    :        (PLBLIN,LBLINT(1))
!     POINTER (PIPTR,IPTR(1))
!
!
      S((-9)) = '-9'
      S((-8)) = '-8'
      S((-7)) = '-7'
      S((-6)) = '-6'
      S((-5)) = '-5'
      S((-4)) = '-4'
      S((-3)) = '-3'
      S((-2)) = '-2'
      S((-1)) = '-1'
      S(0) = '+0'
      S(1) = '+1'
      S(2) = '+2'
      S(3) = '+3'
      S(4) = '+4'
      S(5) = '+5'
      S(6) = '+6'
      S(7) = '+7'
      S(8) = '+8'
      S(9) = '+9'
!
!   Position file at beginning of list of integrals
!
      REWIND (NFILE)
!
!   Initialize
!
      FIRST = .TRUE.
      NMCP = 0
      NINT = 0
!
!   Initial allocation of storage to local arrays
!
      NLABEL = 1
      CALL ALLOC (JLABL, NLABEL, 'JLABL', 'TRSORT')
      CALL ALLOC (XL, NLABEL, 'XL', 'TRSORT')
!
      CALL ALCLLA (IBEG, ILAB, ILAST, ILEFT, IPTCSF, IRIGHT, LBLINT, &
          LLDIM, 1)
      CALL ALCNMA (IPTR, ISLDR, ISLDR1, XSLDR, NMDIM, 1)
!
!   Now the rest of the elements
!
    1 CONTINUE
      READ (NFILE, END=12) IR, IS, NI
      IF (NI > NLABEL) THEN
         CALL RALLOC (JLABL, NI, 'JLABL', 'TRSORT')
         CALL RALLOC (XL, NI, 'XL', 'TRSORT')
         NLABEL = NI
      ENDIF
      READ (NFILE, END=99, ERR=99) (JLABL(I),XL(I),I=1,NI)
      IF (IR==0 .OR. IS==0 .OR. NI==0) GO TO 1
!
      IF (FIRST) THEN
!
!   List is empty
!
         M = 0
         J = 0
!
!   Set up list pointers and insert first element
!
         ICOUNT = 0
    3    CONTINUE
         ICOUNT = ICOUNT + 1
         IF (ICOUNT > NI) GO TO 1
         IF (JLABL(ICOUNT) == 0) GO TO 3
         ILAB(1) = JLABL(ICOUNT)
!
         IRIGHT(1) = 0
         ILEFT(1) = 0
         IBEG(1) = 1
         IPTR(1) = 0
         ILAST(1) = 0
!
         M = 1
         J = 1
!
         FIRST = .FALSE.
!
      ELSE
!
         ICOUNT = 0
!
      ENDIF
!
!   Sort integral list using tree sort
!
!   Take next nonzero element
!
    4 CONTINUE
      ICOUNT = ICOUNT + 1
      IF (ICOUNT > NI) GO TO 1
      JLAB = JLABL(ICOUNT)
      IF (JLAB == 0) GO TO 4
      X = XL(ICOUNT)
!
      M = M + 1
      IF (M > NMDIM) CALL ALCNMA (IPTR, ISLDR, ISLDR1, XSLDR, NMDIM, 2)
      I = 1
!
!   Search for place in tree
!
    5 CONTINUE
      IF (JLAB - ILAB(I) > 0) GO TO 8
      IF (JLAB - ILAB(I) == 0) GO TO 10
      K = IRIGHT(I)
      IF (K /= 0) GO TO 7
      J = J + 1
      IF (J > LLDIM) CALL ALCLLA (IBEG, ILAB, ILAST, ILEFT, IPTCSF, &
                    IRIGHT, LBLINT, LLDIM, 2)
      IRIGHT(I) = J
      GO TO 9
    7 CONTINUE
      I = K
      GO TO 5
    8 CONTINUE
      K = ILEFT(I)
      IF (K /= 0) GO TO 7
      J = J + 1
      IF (J > LLDIM) CALL ALCLLA (IBEG, ILAB, ILAST, ILEFT, IPTCSF, &
                          IRIGHT, LBLINT, LLDIM, 2)
      ILEFT(I) = J
!
!   When found, update list.
!
    9 CONTINUE
      ILAST(J) = I
      IRIGHT(J) = 0
      ILEFT(J) = 0
      IBEG(J) = M
      ILAB(J) = JLAB
      IPTR(M) = 0
      GO TO 4
   10 CONTINUE
      K = IBEG(I)
      L = K
      K = IPTR(L)
      DO WHILE(K /= 0)
         L = K
         K = IPTR(L)
      END DO
      IPTR(L) = M
      IPTR(M) = 0
      GO TO 4
!
!   The end of the CSF-based file has been reached
!
   12 CONTINUE
      IF (.NOT.(FIRST .OR. M==0)) THEN
!
!  Sort is complete. Unpack list
!
         NMCP = M
         NINT = J
         L = 0
         M = 0
         I = 1
!
!   Search for smallest element
!
   13    CONTINUE
         K = IRIGHT(I)
         DO WHILE(K /= 0)
            I = K
            K = IRIGHT(I)
         END DO
!
!   Insert in sorted list
!
   14    CONTINUE
         IF (ILAB(I) == 0) GO TO 16
         L = L + 1
         LBLINT(L) = ILAB(I)
         K = IBEG(I)
!
!   Copy list of pointers to CSF/coefficients into new list
!
         M = M + 1
         ISLDR(M) = K
         K = IPTR(K)
         DO WHILE(K /= 0)
            M = M + 1
            ISLDR(M) = K
            K = IPTR(K)
         END DO
         IPTCSF(L) = M
         ILAB(I) = 0
!
!   Next smallest element is on left of last element
!
         K = ILEFT(I)
         IF (K == 0) GO TO 16
         I = K
         GO TO 13
!
!   If no element on left, next smallest is previous element
!
   16    CONTINUE
         I = ILAST(I)
         IF (I /= 0) GO TO 14
!
!   List is unpacked. Invert CSF/coefficient pointer list to give
!   position list for CSF/coefficients as they are read in
!
         DO I = 1, NMCP
            K = ISLDR(I)
            IPTR(K) = I
         END DO
!
!   Now read CSF pairs and coefficients into correct positions in
!   sorted list
!
         IMCP = 0
         REWIND (NFILE)
   18    CONTINUE
         READ (NFILE, END=20) IR, IS, NI
         IF (IR==0 .OR. IS==0 .OR. NI==0) GO TO 18
         READ (NFILE) (JLABL(I),XL(I),I=1,NI)
         DO I = 1, NI
            IF (JLABL(I) == 0) CYCLE
            IMCP = IMCP + 1
            K = IPTR(IMCP)
            ISLDR(K) = IR
            ISLDR1(K) = IS
!           ISLDR(K) = IR*NCA+IS
            XSLDR(K) = XL(I)
         END DO
         GO TO 18
!
!   The integral-based list is completely known
!
      ENDIF
   20 CONTINUE
      REWIND (NFILE)
!
!  If first set of data open the file and print
!  some data to later be able to identify the file
!
      IF (IBLKI==1 .AND. IBLKF==1) THEN
         J1 = INDEX(NAME(1),' ')
         J2 = INDEX(NAME(2),' ')
         OPEN(UNIT=NFILE2,FILE=NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'//S(&
            KP(JKP))//'T', STATUS='UNKNOWN', FORM='UNFORMATTED', POSITION=&
            'asis')
      ENDIF
      WRITE (NFILE2) IBLKI, IBLKF, NW, NKP
      WRITE (NFILE2) NINT
      IF (NMCP /= 0) THEN
         IF (LPRINT) WRITE (99, 301)
         MLR = 1
         DO I = 1, NINT
            MUP = IPTCSF(I)
            INTS = MUP - MLR + 1
            WRITE (NFILE2) LBLINT(I), INTS
            IF (LPRINT) THEN
               IA = LBLINT(I)/KEY
               IB = MOD(LBLINT(I),KEY)
               WRITE (99, 302) NP(IA), NH(IA), NP(IB), NH(IB)
            ENDIF
            WRITE (NFILE2) (ISLDR(M),ISLDR1(M),XSLDR(M),M=MLR,MUP)
            IF (LPRINT) THEN
               DO M = MLR, MUP
!                  IS = MOD (ISLDR(M),NCA)
                  IS = ISLDR1(M)
!                 IR = ISLDR(M)/NCA
                  IR = ISLDR(M)
                  WRITE (99, 303) IR, IS, XSLDR(M)
               END DO
            ENDIF
            MLR = MUP + 1
         END DO
!
      ENDIF
!
!   Deallocate storage for local arrays
!
      CALL DALLOC (JLABL, 'JLABJ', 'TRSORT')
      CALL DALLOC (XL, 'XL', 'TRSORT')
      CALL ALCLLA (IBEG, ILAB, ILAST, ILEFT, IPTCSF, IRIGHT, LBLINT, &
                    LLDIM, 3)
      CALL ALCNMA (IPTR, ISLDR, ISLDR1, XSLDR, NMDIM, 3)
!
      RETURN
!
!   Error handling
!
   99 CONTINUE
      WRITE (6, *) 'TRSORT: Error reading CSF-based file.'
      STOP
!
  301 FORMAT(/,/,'  k'/,' d  (rs)  Coefficients:'/,'  ab'/,/,/,&
         '    a     b          r         s  Coefficient'/)
  302 FORMAT(2(2X,I2,A2))
  303 FORMAT(14X,1I6,2X,1I6,2X,1P,1D22.15)
      RETURN
!
      END SUBROUTINE TRSORT
