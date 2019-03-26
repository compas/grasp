!***********************************************************************
!                                                                      *
      SUBROUTINE ORBORD(N)
!                                                                      *
!   THIS ROUTINE DOES NOTHING!

!   This subroutine checks the orbital order. Normal and reversed      *
!   order can be used. If reveresed order JA and JB must be            *
!   permuted in the subroutine ti1tv                                   *
!                                                                      *
!   Written by Per Jonsson                Last revision: Feb    1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE orb_C
      USE biorb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------

      RETURN
      END SUBROUTINE ORBORD
