!***********************************************************************
!                                                                      *
      SUBROUTINE FACTT
!                                                                      *
!   Calculates the logs  of factorials required by the Racah coeffi-   *
!   cient routine DRACAH.                                              *
!                                                                      *
!   Written by N S Scott                    Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:26   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE FACTS_C
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: X
!-----------------------------------------------
!
      GAM(1) = 1.0D00
      GAM(2) = 1.0D00
      X = 2.0D00
!
      DO I = 3, 30
         GAM(I) = GAM(I-1)*X
         X = X + 1.0D00
      END DO
!
      DO I = 1, 30
         GAM(I) = LOG(GAM(I))
      END DO
!
      X = 3.0D01
!
      DO I = 31, MFACT
         GAM(I) = GAM(I-1) + LOG(X)
         X = X + 1.0D00
      END DO
!
      RETURN
      END SUBROUTINE FACTT
