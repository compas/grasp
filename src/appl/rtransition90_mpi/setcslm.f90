!***********************************************************************
!                                                                      *
      SUBROUTINE SETCSLM 
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to: [RCI92]: LENGTH, LODCSL, OPENFL.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      USE lodcslm_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, IOS, NCORE 
      LOGICAL :: FOUND 
      CHARACTER :: FILNAM*256, RECORD*15, DEFNAM*11, FORM*11, STATUS*3 
!-----------------------------------------------
!
!
!   The  .csl  file is FORMATTED; it must exist
!
      DEFNAM = 'SLASK' 
      FORM = 'FORMATTED' 
      STATUS = 'OLD' 
!
!   Look for  grasp92.csl
!
      INQUIRE(FILE=DEFNAM, EXIST=FOUND) 
!
      IF (FOUND) THEN 
         FILNAM = DEFNAM 
      ELSE 
         WRITE (6, *) 'rcsl.inp does not exist' 
         STOP  
      ENDIF 
!
      CALL OPENFL (21, FILNAM, FORM, STATUS, IERR) 
      IF (IERR == 1) THEN 
         WRITE (6, *) 'Error when opening rcsl.inp' 
         STOP  
      ENDIF 
!
!   Check the first record of the file; if not as expected, try again
!
      READ (21, '(1A15)', IOSTAT=IOS) RECORD 
      IF (IOS/=0 .OR. RECORD(1:15)/='Core subshells:') THEN 
         WRITE (6, *) 'Not a Configuration Symmetry List File;' 
         CLOSE(21, STATUS='DELETE') 
      ENDIF 
!
!   Load data from the  .csl  file
!
      CALL LODCSLM (NCORE) 
!
!   Close the  .csl  file
!
      CLOSE(21, STATUS='DELETE') 
!
      RETURN  
      END SUBROUTINE SETCSLM 
