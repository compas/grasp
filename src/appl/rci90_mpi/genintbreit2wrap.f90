!***********************************************************************
!                                                                      *
        SUBROUTINE genintbreit2wrap (myid, nprocs, j2max)

!   Written by Per Jonsson                Last revision: October 2014  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE bilst_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE genintbreit2_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: myid, nprocs
      INTEGER             :: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: N
!-----------------------------------------------
!
      CALL genintbreit2 ((myid), (nprocs), N, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      CALL gisummpi (INDTP2, N)
      CALL gdsummpi (VALTP2, N)

      RETURN
      END SUBROUTINE genintbreit2wrap
