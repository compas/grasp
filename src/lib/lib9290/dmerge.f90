!**************************************************************************
!
      SUBROUTINE DMERGE(N, DB, DC, IDY, DA, DCONST, DL)
!-----------------------------------------------
!
!  this merge version has the advantage of loading da(i)
!  and idy(i) only once.
!
!**************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:20   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)    :: N
      REAL(DOUBLE), INTENT(IN)    :: DCONST
      REAL(DOUBLE), INTENT(OUT)   :: DL
      INTEGER, DIMENSION(N), INTENT(IN) :: IDY
      REAL(DOUBLE), DIMENSION(*), INTENT(IN) :: DB
      REAL(DOUBLE), DIMENSION(*), INTENT(INOUT) :: DC
      REAL(DOUBLE), DIMENSION(N), INTENT(IN) :: DA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: DSUM
!-----------------------------------------------

      DSUM = 0.0
      DO I = 1, N
         DSUM = DSUM + DA(I)*DB(IDY(I))
         DC(IDY(I)) = DC(IDY(I)) + DCONST*DA(I)
      END DO
      DL = DSUM

      RETURN
      END SUBROUTINE DMERGE
