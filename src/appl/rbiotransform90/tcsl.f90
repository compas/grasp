!***********************************************************************
!                                                                      *
      SUBROUTINE TCSL(N)
!                                                                      *
!   This subroutine transfers data to the initial and final state      *
!   common blocks                                                      *
!                                                                      *
!   Written by Per Jonsson                Last revision: June   1996   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE biorb_C
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------
!
!  Initial state commons
!
!
!  Final state commons
!
!
      IF (N == 1) THEN
         NCFII = NCF
         NWII = NW
         NPII(:NW) = NP(:NW)
         NAKII(:NW) = NAK(:NW)
         NKLII(:NW) = NKL(:NW)
         NKJII(:NW) = NKJ(:NW)
         NHII(:NW) = NH(:NW)
      ELSE
         NCFFF = NCF
         NWFF = NW
         NPFF(:NW) = NP(:NW)
         NAKFF(:NW) = NAK(:NW)
         NKLFF(:NW) = NKL(:NW)
         NKJFF(:NW) = NKJ(:NW)
         NHFF(:NW) = NH(:NW)
      ENDIF

      RETURN
      END SUBROUTINE TCSL
