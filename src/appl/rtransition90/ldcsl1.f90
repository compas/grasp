!***********************************************************************
!                                                                      *
      SUBROUTINE LDCSL1 (NCORER,NAME)
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
!                                                                      *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE memory_man
      USE def_C
      USE biorb_C
      USE orb_C
      USE stat_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      USE lodcsl_I
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ncorer
      CHARACTER(LEN=24) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: i, j, ios, ncore
       CHARACTER(LEN=15) :: record
!
!   Check the first record of the file; if not as expected, try again
!
      J = INDEX(NAME,' ')
      OPEN (UNIT = 21,FILE=NAME(1:J-1)//'.c',FORM='FORMATTED',          &
     &     STATUS='OLD')

      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.                                             &
     &    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (21)
         STOP
      ENDIF
!
!   Load data from the  .csl  file
!
      CALL LODCSL (NCORE)
!
!   Close the  .csl  file
!
      CLOSE (21)
!
!   Check if the core should be redefined
!
      NCORER = NCORE
!      DO 3 I = NCORE+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORER = NCORER+1
!    3 CONTINUE
!
!   Copy data to reference arrays
!
    4 NELECR = NELEC
      NELECII = NELEC
      NCFR = NCF
      NWR = NW
!
      CALL ALLOC (IQAR, NNNW, NCF, 'IQAR', 'LDCSL1' )
      CALL ALLOC (JQSAR,NNNW,3,NCF, 'JQSAR', 'LDCSL1')
      CALL ALLOC (JCUPAR,NNNW, NCF, 'JCUPAR', 'LDCSL1')
!
      DO 6 I = 1,NCF
         DO 5 J = 1,NNNW
            IQAR(J,I) = IQA(J,I)
            JQSAR(J,1,I) = JQSA(J,1,I)
            JQSAR(J,2,I) = JQSA(J,2,I)
            JQSAR(J,3,I) = JQSA(J,3,I)
            JCUPAR(J,I) = JCUPA(J,I)
    5    CONTINUE
    6 CONTINUE
!
      NWII=NW
      NCFII=NCF
      DO 7 I = 1,NW
         NHII(I) = NH(I)
         NPII(I) = NP(I)
         NAKII(I) = NAK(I)
         NKLII(I) = NKL(I)
         NKJII(I) = NKJ(I)
         NPR(I) = NP(I)
         NAKR(I) = NAK(I)
         NHR(I) = NH(I)
    7 CONTINUE
!
!   Deallocate storage so that LODCSL can restart correctly
!
      CALL DALLOC (IQA, 'IQA', 'LDCSL1')
      CALL DALLOC (JQSA, 'JQSA', 'LDCSL1')
      CALL DALLOC (JCUPA, 'JCUPA', 'LDCSL1')
!
      CLOSE (21)

      RETURN
      END
