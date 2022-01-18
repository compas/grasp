!***********************************************************************
!                                                                      *
      SUBROUTINE SETXV(J)
!                                                                      *
!   This  subprogram  sets up the inhomogeneous terms for the varia-   *
!   tion equations.                                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE grid_C
      USE int_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: HHC
!-----------------------------------------------
!
      HHC = 0.5D00*H/C
!
!   Set up arrays TF and TG
!
      DO I = 1, N
         XU(I) = -QF(I,J)*HHC*RP(I)
         XV(I) = PF(I,J)*HHC*RP(I)
      END DO
!
      XU(:N-1) = XU(2:N) + XU(:N-1)
      XV(:N-1) = XV(2:N) + XV(:N-1)
!
      RETURN
!
      END SUBROUTINE SETXV
