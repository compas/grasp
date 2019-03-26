!***********************************************************************
!                                                                      *
      SUBROUTINE SETXZ(J)
!                                                                      *
!   This subprogram sets the inhomogeneous terms to zero.              *
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
      USE grid_C
      USE int_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
!-----------------------------------------------
!
      XU(:N) = 0.0D00
      XV(:N) = 0.0D00
!
      RETURN
      END SUBROUTINE SETXZ
