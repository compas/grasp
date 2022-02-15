!***********************************************************************
!                                                                      *
      SUBROUTINE MAXARR(J)
!                                                                      *
!   This subroutine finds the least self-consistent orbital            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
! J initialized to zero
! XHH 1997.02.14
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE fixd_C
      USE orb_C
      USE scf_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: DLRGST
!-----------------------------------------------
!
!
      J = 0
      DLRGST = 0.D0
      DO I = 1, NW
         IF (LFIX(I)) CYCLE
         IF (SCNSTY(I) <= DLRGST) CYCLE
         DLRGST = SCNSTY(I)
         J = I
      END DO
!
      RETURN
      END SUBROUTINE MAXARR
