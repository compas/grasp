!***********************************************************************
!
      SUBROUTINE SETCSH(NFILE, NAME, NCORE) 
!
!   Open, check the CSL file and load the load (via lodcsh) data from
!   the header lines. It is designed to replace all kinds of "setcsl"
!   routines within GRASP packages.
!
!   Routines called: lodcsh
!
!   ncore is output parameter
!
!   Written by Xinghone He                Last revision: 23 Dec 1997
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:30   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE IOUNIT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lodcsh_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NFILE 
      INTEGER  :: NCORE 
      CHARACTER (LEN = *), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENGTH, IOS 
      LOGICAL :: FOUND 
      CHARACTER :: RECORD*15 
!-----------------------------------------------
!
!     ...Locals
 
      LENGTH = LEN_TRIM(NAME) 
 
      INQUIRE(FILE=NAME, EXIST=FOUND) 
      IF (.NOT.FOUND) THEN 
         WRITE (ISTDE, *) NAME(1:LENGTH), ' does not exist' 
         STOP  
      ENDIF 
 
      INQUIRE(UNIT=NFILE, OPENED=FOUND) 
      IF (FOUND) THEN 
         WRITE (ISTDE, *) 'Unit ', NFILE, ' has been used elsewhere' 
         STOP  
      ENDIF 
 
      OPEN(NFILE, FILE=NAME, STATUS='OLD', IOSTAT=IOS, POSITION='asis') 
      IF (IOS /= 0) THEN 
         WRITE (ISTDE, *) 'Error when opening ', NAME 
         STOP  
      ENDIF 
!
!   Check the first record of the file; if not as expected, try again
!
      READ (NFILE, '(1A15)', IOSTAT=IOS) RECORD 
 
      IF (IOS/=0 .OR. RECORD(1:15)/='Core subshells:') THEN 
         WRITE (ISTDE, *) 'Not a Configuration Symmetry List File;' 
         CLOSE(NFILE) 
         STOP  
      ENDIF 
!
!   Load data from the  .csl  file
!
      CALL LODCSH (NFILE, NCORE) 
 
      RETURN  
      END SUBROUTINE SETCSH 
