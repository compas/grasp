!***********************************************************************

      SUBROUTINE SETSUM(NAME)
!
!   Open the  .csum  file on stream 24.
!   Xinghong He                                          10 Jun 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) , INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IERR
!-----------------------------------------------------------------------

      K = INDEX(NAME,' ')

      CALL OPENFL (24, NAME(1:K-1)//'.csum', 'FORMATTED', 'UNKNOWN', IERR)

      RETURN
      END SUBROUTINE SETSUM
