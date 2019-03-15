!***********************************************************************
!                                                                      *
      SUBROUTINE LDCSL2(NCORE, NAME)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE orb_C
      USE biorb_C

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lodcsl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: NCORE
      CHARACTER , INTENT(IN) :: NAME*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, IOS, NCORER, I
      CHARACTER :: RECORD*15
!-----------------------------------------------
!
!
!     Common relevant for the radial wave functions of the final state.
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      J = INDEX(NAME,' ')
      OPEN(UNIT=21, FILE=NAME(1:J-1)//'.c', FORM='FORMATTED', STATUS='OLD', &
         POSITION='asis')
!
!   Check the first record of the file; if not as expected, try again
!
      READ (21, '(1A15)', IOSTAT=IOS) RECORD
      IF (IOS/=0 .OR. RECORD(1:15)/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         CLOSE(21)
         STOP
      ENDIF
!
!   Load data from the  .csl  file
!
      CALL LODCSL (NCORER)
!
!   Close the  .csl  file
!
      CLOSE(21)
!
!   Check if the core should be redefined
!
      NCORE = NCORER
!      DO 3 I = NCORER+1,NW
!         IFULLI = NKJ(I)+1
!         DO 2 J = 1,NCF
!            IF (IQ (I,J) .NE. IFULLI) GOTO 4
!    2    CONTINUE
!         CALL CONVRT (NP(I),RECORD,LENTH)
!         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
!         PRINT *, ' in all CSFs; including this'
!         PRINT *, ' subshell in the core;'
!         NCORE = NCORE+1
!    3 CONTINUE
!
      NELECFF = NELEC
      NWFF = NW
      NCFFF = NCF
      NHFF(:NW) = NH(:NW)
      NPFF(:NW) = NP(:NW)
      NAKFF(:NW) = NAK(:NW)
      NKLFF(:NW) = NKL(:NW)
      NKJFF(:NW) = NKJ(:NW)
      RETURN
      END SUBROUTINE LDCSL2
