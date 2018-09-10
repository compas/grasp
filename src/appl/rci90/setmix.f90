!***********************************************************************
!                                                                      *
      SUBROUTINE SETMIX(NAME, IDBLK) 
!                                                                      *
!   Opens the  .mix  file on stream 25; writes a header to this file;  *
!   calls LODMIX to interactively determine the eigenpairs required.   *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                                       *
!               [RCI92]: LODMIX.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Modified by Xinghong He               Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      USE lodmix_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) , INTENT(IN) :: NAME 
      CHARACTER(LEN=8), DIMENSION(*) :: IDBLK 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*11, PARAMETER :: FORM = 'UNFORMATTED' 
      CHARACTER*7, PARAMETER :: STATUS = 'UNKNOWN' 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IERR 
!-----------------------------------------------
 
      K = INDEX(NAME,' ') 
      CALL OPENFL (25, NAME(1:K-1)//'.cm', FORM, STATUS, IERR) 
      IF (IERR /= 0) STOP 'setmix: Error when opening .cm file' 
 
      WRITE (25) 'G92MIX' 
 
      CALL LODMIX (IDBLK) 
 
      RETURN  
      END SUBROUTINE SETMIX 
