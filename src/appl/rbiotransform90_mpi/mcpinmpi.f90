!***********************************************************************
!                                                                      *
      SUBROUTINE MCPIN (NAME,startdir,IK,NTESTG,INPCI)
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
!   Bug corrected 2005-10-18                                           *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB, NNNW
      USE memory_man
      USE def_C
      USE eigv_C
      USE foparm_C
      USE mcp_C
      USE orb_C,           ONLY: ncf, nw
      USE prnt_C
      USE syma_C
      USE stat_C
      USE jqjc_C
      USE orbord_C
      USE mcpdata_C
      USE sbdat_C,         ONLY: KAMAX,NSHLII,NSHLFF,NSHLPPII,NSHLPPFF,&
                                 NINII,NINFF,NSHLPII,NSHLPFF,NAKINVII, &
                                 NAKINVFF,NLMAX
      USE sbdat1_C,        ONLY: nshlp, nshlpp
      USE cimat_C,         ONLY: cici, cfci
      USE blk_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE tnsrjj_I
!      USE cor_I
!      USE cord_I
!      USE qqsort_I
!      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=24), INTENT(IN)  :: name
      CHARACTER(LEN=128), INTENT(IN) :: startdir
      INTEGER, INTENT(IN) :: ik, ntestg, inpci
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!   Set the encoding key and its square
      INTEGER, PARAMETER :: KEY = KEYORB, KEYSQ = KEY*KEY, nvmax=100
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
      INTEGER, PARAMETER :: NF = 200

      LOGICAL ::  F0INT,LINCR,RESTRT,COMP,AVAIL
      CHARACTER(LEN=20) :: CNUM
      CHARACTER(LEN=2)  :: CK
!
      REAL(DOUBLE), DIMENSION(NNNW) :: tshell
      REAL(DOUBLE), DIMENSION(NVMAX) ::SC
      INTEGER, DIMENSION(NNNW) :: nakinv
      INTEGER, DIMENSION(NLMAX) :: llistt, nshl, ninl
      REAL(DOUBLE), DIMENSION(20*NLMAX*NLMAX) :: cirot
      REAL(DOUBLE), DIMENSION(10000) :: EVSC
!
      REAL(DOUBLE), DIMENSION(:), pointer :: scr, ciout
      INTEGER :: ntest, ntestL, i, j, iblk, ncfd, nwd, kamaxd,  &
                 l, iioff, nvecsize
!
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 00

!
!   Set up data for the initial and final state case respectively
!
      IF (IK.EQ.1) THEN
        DO I = 1,NLMAX
          NSHL(I) = NSHLII(I)
          NINL(I) = NINII(I)
        ENDDO

        DO I = 1,NNNW
          NAKINV(I) = NAKINVII(I)
        ENDDO

        DO J = 1,NLMAX
          DO I = 1,NLMAX
            NSHLP(I,J) = NSHLPII(I,J)
          ENDDO
        ENDDO

        DO J = 1,NNNW
          DO I = 1,NLMAX
            NSHLPP(I,J) = NSHLPPII(I,J)
          ENDDO
        ENDDO

        DO I = 1,20*NLMAX*NLMAX
          CIROT(I) = CICI(I)
        ENDDO

      ELSEIF (IK.EQ.2) THEN

        DO I = 1,NLMAX
          NSHL(I) = NSHLFF(I)
          NINL(I) = NINFF(I)
        ENDDO

        DO I = 1,NNNW
          NAKINV(I) = NAKINVFF(I)
        ENDDO

        DO J = 1,NLMAX
          DO I = 1,NLMAX
            NSHLP(I,J) = NSHLPFF(I,J)
          ENDDO
        ENDDO

        DO J = 1,NNNW
          DO I = 1,NLMAX
            NSHLPP(I,J) = NSHLPPFF(I,J)
          ENDDO
        ENDDO

        DO I = 1,20*NLMAX*NLMAX
          CIROT(I) = CFCI(I)
        ENDDO
      ENDIF

      REWIND (NF)
      DO 1000 IBLK = 1, NBLOCK
!
!   Read the CI vectors
!
      CALL GETMIX(trim(startdir)//'/'//NAME,INPCI,IBLK)
!
!   Allocate memory for use in CITRAG
!
      CALL alloc (scr, ncf*nvec, 'scr', 'MCPIN')
      CALL alloc (ciout, ncf*nvec, 'CIOUT', 'MCPIN')
!
!   if data available on file read the datafile. Else sort the
!   MCP data into inegral based lists.
!
      READ(NF) NCFD,NWD,KAMAXD
      DO L = 1,KAMAX

!************
!
!. Offset for given L in shell matrices
!
          IF (L.EQ.1) THEN
            IIOFF = 1
          ELSE
            IIOFF = IIOFF + NSHL(L-1)** 2
          END IF
!  Corrected PER J

!**************

          READ(NF) NINTG,NCOEFF
          IF(NCOEFF*NINTG.EQ.0) CYCLE
!
!   Allocate memory. Note that in this case we can,
!   since the data is sorted, allocate less
!   memory than is done in qqsort.
!
          CALL ALLOC (JANN,NCOEFF, 'JANN', 'MCPIN')
          CALL ALLOC (JBNN,NCOEFF, 'JBNN', 'MCPIN')
          CALL ALLOC (INTGRL,NINTG, 'INTGRL', 'MCPIN')
          CALL ALLOC (CNN,NCOEFF, 'CNN', 'MCPIN')
          CALL ALLOC (INTPTR,NINTG, 'INTPTR', 'MCPIN')

          DO I = 1,NINTG
            READ(NF) INTGRL(I),INTPTR(I)
          ENDDO
          DO I = 1,NCOEFF
            READ(NF) CNN(I),JANN(I),JBNN(I)
          ENDDO
!*
!*. Offset for given L in shell matrices
!*
!        IF (L.EQ.1) THEN
!          IIOFF = 1
!        ELSE
!          IIOFF = IIOFF + NSHL(L-1)** 2
!        END IF
!   BUG THE OFFSET SHOULD BE BEFORE THE CYCLE STATEMENT!!
!
! Transform for given L
!
        IF (NSHL(L).GT.0)                                             &
!        CALL CITRAG(EVEC,NCF,NVEC,L,NSHL(L),                   &
        CALL CITRAG(EVEC(1:500),NCF,NVEC,L,NSHL(L),                   &
                    CIROT(IIOFF),NINL(L),NTESTL,CIOUT, SCR )
!
!   Deallocate storage for this kappa.
!
        CALL DALLOC (JANN, 'JANN', 'MCPIN')
        CALL DALLOC (JBNN, 'JBNN', 'MCPIN')
        CALL DALLOC (INTGRL, 'INTGRL', 'MCPIN')
        CALL DALLOC (CNN, 'CNN', 'MCPIN')
        CALL DALLOC (INTPTR, 'INTPRT', 'MCPIN')

      ENDDO

      CALL dalloc (scr, 'SCR', 'MCPIN')
      CALL dalloc (ciout, 'CIOUT', 'MCPIN')
!
!   Write the rotated CI vectors on file
!
      IF(myid .eq. 0) THEN
      IF(IBLK.EQ.1) THEN
        J = INDEX(NAME,' ')
        IF (INPCI.EQ.0) THEN
          OPEN (UNIT=31,FILE=trim(startdir)//'/'//NAME(1:J-1)//'.cbm', &
                FORM='UNFORMATTED',STATUS='UNKNOWN')
        ELSE
          OPEN (UNIT=31,FILE=trim(startdir)//'/'//NAME(1:J-1)//'.bm',  &
                FORM='UNFORMATTED',STATUS='UNKNOWN')
        ENDIF

        WRITE(31) 'G92MIX'
        WRITE(31) NELEC,NCFTOT,NW,NVECTOT,NVECSIZE,NBLOCK
      ENDIF
         WRITE(31) IBLK,NCF,NVEC,IATJPO(1),IASPAR(1)
         WRITE(31) (IVEC(I),I = 1,NVEC)
         WRITE(31) EAV,(EVAL(I),I = 1,NVEC)
         WRITE(31) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
      endif !myid=0

      CALL DALLOC (EVAL, 'EVAL', 'MCPIN')
      CALL DALLOC (EVEC, 'EVEC', 'MCPIN')
      CALL DALLOC (IVEC, 'IVEC', 'MCPIN')
      CALL DALLOC (IATJPO, 'IATJPO', 'MCPIN')
      CALL DALLOC (IASPAR, 'IASPAR', 'MCPIN')
 1000 CONTINUE
      CLOSE (NF)
!
!   Close the  mixing  file
!
      CLOSE (30)
      CLOSE (31)
!
      RETURN
      END SUBROUTINE MCPIN
