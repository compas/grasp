!***********************************************************************
!                                                                      *
      SUBROUTINE SETTMP(NB, KMAX, FILEHEAD)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      USE openfl_I
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE iounit_C
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NB
      INTEGER , INTENT(IN) :: KMAX
      CHARACTER , INTENT(IN) :: FILEHEAD*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, LCK, IERR, I, LNG
      CHARACTER :: CK*2
!-----------------------------------------------
!
! All files  filehead.XX  are UNFORMATTED;
!
      LNG = LEN_TRIM(FILEHEAD)
      DO K = 30, 32 + KMAX
         CALL CONVRT (K, CK, LCK)
         CALL OPENFL (K, FILEHEAD(1:LNG)//'.'//CK(1:2), 'UNFORMATTED', &
            'UNKNOWN', IERR)
         IF (IERR == 0) CYCLE
         DO I = 30, K
            CLOSE(I)
         END DO
         WRITE (ISTDE, *) 'Error when opening the tmp files'
         STOP
      END DO

      DO K = 30, 32 + KMAX
         WRITE (K) 'MCP', NB
      END DO

      RETURN
      END SUBROUTINE SETTMP
