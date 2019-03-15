!***********************************************************************
!                                                                      *
      SUBROUTINE MCTOUT(IOPAR, JKP, NAME)
!                                                                      *
!   This routine loads  coefficients with parity and  rank specified   *
!   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
!   (+/- 1) and is determined from the sign of  KP(JKP).               *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:29:28   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB, NNNW
      USE memory_man
      USE blk_C
      USE debug_C,         ONLY: LDBPA
      USE decide_C
      USE foparm_C
      USE orb_C
      USE OFFD_C
      USE OSC_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE angdata_I
      USE itjpo_I
      USE oneparticlejj_I
      USE trsort_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IOPAR
      INTEGER  :: JKP
      CHARACTER  :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
      INTEGER, PARAMETER :: NFILE = 93
      INTEGER, PARAMETER :: NFILE1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NFILE2, IBLKI, IBLKF, NMCT, NLABEL, NCFI0, NCFF0, IC, IR, NCR&
         , IA, IB, NEWSIZ, I
      LOGICAL :: AVAIL
      INTEGER, DIMENSION(:), pointer :: label
      REAL(DOUBLE), DIMENSION(:), pointer:: coeff
      REAL(DOUBLE), DIMENSION(NNNW) :: tshell
!-----------------------------------------------
!
!
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
!
!   Check if angular data is available on file
!
      NFILE2 = NFILE1 + JKP
      CALL ANGDATA (NAME, AVAIL, JKP, NFILE2)

!
!   If angular data is not available open the scratch file to store the
!   coefficients; position file
!   to beginning
!
!
      LK = ABS(KP(JKP))
      IOPAR = ISIGN(1,KP(JKP))
      IF (AVAIL) RETURN
      WRITE (6, *) 'LK,IOPAR,from MCTOUT'
      WRITE (6, *) LK, IOPAR
!
! Start of the block loops
      DO IBLKI = 1, NBLOCKI
         DO IBLKF = 1, NBLOCKF
!        OPEN (NFILE,STATUS = 'new', FORM = 'UNFORMATTED')
!        Sometimes, when there has been an error, status need to be "unknown"
            OPEN(NFILE,STATUS='unknown',FORM='UNFORMATTED')
            NMCT = 0
!
!
!   If angular data is not available
!   Generate MCT coefficients for given rank/parity combination and
!   store them by CSF on NFILE
!
!
!   Allocate storage to buffer arrays
!
            NLABEL = 32
            CALL ALLOC (LABEL, NLABEL, 'LABEL', 'MCTOUT')
            CALL ALLOC (COEFF, NLABEL, 'COEFF', 'MCTOUT')

            IF (IBLKI == 1) THEN
               NCFI0 = 1
            ELSE
               NCFI0 = NCFI(IBLKI-1) + 1
            ENDIF
!
            IF (IBLKF == 1) THEN
               NCFF0 = NCFI(NBLOCKI) + 1
            ELSE
               NCFF0 = NCFI(NBLOCKI) + NCFF(IBLKF-1) + 1
            ENDIF
            DO IC = NCFI0, NCFI(IBLKI)
               DO IR = NCFF0, NCFI(NBLOCKI) + NCFF(IBLKF)
!
!           IR = IC
!
                  NCR = 0

!
!   In many case one is interested only in M1 and E2 transitions between
!   levels with different J values. If this is the case then the do check
!   on the J quantum numbers of the CSFs before calling TNSRJJ.
!
                  IF (KP(JKP)==1 .AND. NOFFD1==1) THEN
                     IF (ITJPO(IC) == ITJPO(IR)) CYCLE
                  ENDIF
                  IF (KP(JKP)==2 .AND. NOFFD2==1) THEN
                     IF (ITJPO(IC) == ITJPO(IR)) CYCLE
                  ENDIF
!          if(ispar(ic)*ispar(ir)*iopar.ne.1.
!    &        or.itrig(itjpo(ic),itjpo(ir),2*lk+1).ne.1) go to 13
!          if(ichkq1(IC,IR).eq.0) go to 13
                  CALL ONEPARTICLEJJ(LK,IOPAR,IC,IR,IA,IB,TSHELL)
                  IF (IA /= 0) THEN
                     IF (IA == IB) THEN
                        DO IA = 1, NW
                           IF (ABS(TSHELL(IA)) <= CUTOFF) CYCLE
                           NCR = NCR + 1
                           IF (NCR > NLABEL) THEN
                              NEWSIZ = 2*NLABEL
                              CALL RALLOC (LABEL,  NEWSIZ, 'LABEL', 'MCTOUT')
                              CALL RALLOC (COEFF,  NEWSIZ, 'COEFF', 'MCTOUT')
                              NLABEL = NEWSIZ
                           ENDIF
                           LABEL(NCR) = IA*KEYORB + IA
                           COEFF(NCR) = TSHELL(IA)
                        END DO
                     ELSE
                        IF (ABS(TSHELL(1)) > CUTOFF) THEN
                           NCR = NCR + 1
                           IF (NCR > NLABEL) THEN
                              NEWSIZ = 2*NLABEL
                              CALL RALLOC (LABEL,  NEWSIZ, 'LABEL', 'MCTOUT')
                              CALL RALLOC (COEFF,  NEWSIZ, 'COEFF', 'MCTOUT')
                              NLABEL = NEWSIZ
                           ENDIF
                           LABEL(NCR) = IA*KEYORB + IB
                           COEFF(NCR) = TSHELL(1)
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (NCR <= 0) CYCLE
                  WRITE (NFILE) IC - NCFI0 + 1, IR - NCFF0 + 1, NCR
!             WRITE (NFILE) IC,IR,NCR
                  WRITE (NFILE) (LABEL(I),COEFF(I),I=1,NCR)
                  NMCT = NMCT + NCR

!
               END DO
            END DO
!
!   Deallocate storage for buffer arrays
!
            CALL DALLOC (LABEL,  'LABEL', 'MCTOUT')
            CALL DALLOC (COEFF,  'COEFF', 'MCTOUT')
!
            WRITE (*, 301) NMCT, LK, IOPAR
!
!   Sort the MCT coefficients by integral labels
!
            CALL TRSORT (NAME, NFILE, NFILE2, LDBPA(2), JKP, IBLKI, IBLKF)
            CLOSE(NFILE, STATUS='delete')
! end of the loops for blocks
         END DO
      END DO
!
!   Read the data back as required by OSCL conventions
!
      REWIND (NFILE2)
      RETURN
!
  301 FORMAT(/,/,/,1X,I8,' MCT coefficients generated for rank ',I2,&
         ' and parity ',I2,/,/)
      RETURN
!
      END SUBROUTINE MCTOUT
