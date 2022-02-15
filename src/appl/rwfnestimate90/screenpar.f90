

!***********************************************************************
      SUBROUTINE SCREENPAR(NCORE)
!
! Purpose:
!   Compute hydrogenic screen parameters
!
! Input:
!   ncore
!   in the common - nw, nkj()
!
! Output:
!   in the common - sigma()
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE ORB_C
      USE HYDPAR_C, ONLY: SIGMA
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NCORE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NELECTRON, I
!-----------------------------------------------

       !...Core orbitals
      NELECTRON = 0
      DO I = 1, NCORE
         SIGMA(I) = NELECTRON + (NKJ(I)+1)/2
         NELECTRON = NELECTRON + NKJ(I) + 1
      END DO

       !...Peel orbitals
      SIGMA(NCORE+1:NW) = NELECTRON

      RETURN
      END SUBROUTINE SCREENPAR
