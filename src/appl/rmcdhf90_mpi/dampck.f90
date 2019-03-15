!***********************************************************************
!                                                                      *
      SUBROUTINE DAMPCK(IPR, J, ED1, ED2)
!                                                                      *
!   This subroutine determines the damping factor appropriate to the   *
!   present  orbital. The algorithm is taken from C Froese Fischer's   *
!   program MCHF.                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!   Modified by C. Froese Fischer           Last update: 07 Apr 2009   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE damp_C
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: IPR
      INTEGER, INTENT(IN) :: J
      REAL(DOUBLE), INTENT(INOUT) :: ED1
      REAL(DOUBLE), INTENT(INOUT) :: ED2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: ADAPTV
!-----------------------------------------------
!
!   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
!   if and only if ODAMP(J) .GE. 0.0
!
      ADAPTV = ODAMP(J) >= 0.0D00
!
      IF (IPR /= J) THEN
         IF (ADAPTV) ODAMP(J) = 0.75D00*ODAMP(J)
      ELSE
         ED2 = ED2-E(J)
         IF (ADAPTV) THEN
            IF (ED1*ED2 > 0.0D00) THEN
               ODAMP(J) = 0.75D00*ODAMP(J)
            ELSE
               ODAMP(J) = 0.25D00+0.75D00*ODAMP(J)
            ENDIF
         ENDIF
      ENDIF
!
      IF (IPR == J) THEN
         ED1 = ED2
      ELSE
         ED1 = ED2-E(J)
      ENDIF
!
      IPR = J
      RETURN
      END SUBROUTINE DAMPCK
