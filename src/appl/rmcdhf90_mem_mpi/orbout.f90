!***********************************************************************
      SUBROUTINE ORBOUT(RWFFILE2)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE grid_C
      USE orb_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER , INTENT(IN) :: RWFFILE2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, MFJ, I
!-----------------------------------------------
!
      OPEN(23, FILE=RWFFILE2, STATUS='UNKNOWN', FORM='UNFORMATTED', POSITION=&
         'asis')
      WRITE (23) 'G92RWF'
      DO J = 1, NW
         MFJ = MF(J)
         WRITE (23) NP(J), NAK(J), E(J), MFJ
         WRITE (23) PZ(J), (PF(I,J),I=1,MFJ), (QF(I,J),I=1,MFJ)
         WRITE (23) (R(I),I=1,MFJ)               ! This is a waste of resources
      END DO
      CLOSE(23)

      RETURN
      END SUBROUTINE ORBOUT
