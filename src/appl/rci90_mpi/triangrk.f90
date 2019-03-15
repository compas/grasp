!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TRIANGRK (LA, K, LB)
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: LA
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: LB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!   Perform the triangularity check
!
      IF (MOD(K + LA + LB,2) /= 0) THEN
         TRIANGRK = .FALSE.
      ELSE
         IF (ABS(LA - LB) > K) THEN
            TRIANGRK = .FALSE.
         ELSE IF (LA + LB < K) THEN
            TRIANGRK = .FALSE.
         ELSE
            TRIANGRK = .TRUE.
         ENDIF
      ENDIF

      RETURN
      END FUNCTION TRIANGRK
