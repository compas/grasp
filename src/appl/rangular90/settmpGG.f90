!***********************************************************************
!                                                                      *
      SUBROUTINE SETTMPGG(NB, K, FILEHEAD) 
!***********************************************************************
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
!***********************************************************************
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
      INTEGER, INTENT(IN) :: NB 
      INTEGER, INTENT(IN) :: K 
      CHARACTER, INTENT(IN) :: FILEHEAD*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: LCK, IERR, LNG 
      CHARACTER :: CK*2 
!-----------------------------------------------
!
! All files  filehead.XX  are UNFORMATTED;
!
      LNG = LEN_TRIM(FILEHEAD)
      CALL CONVRT (K, CK, LCK)
      CALL OPENFL (K, FILEHEAD(1:LNG)//'.'//CK(1:2), 'UNFORMATTED', &
            'UNKNOWN', IERR)
      IF (IERR /= 0) THEN
         CLOSE(K)
         WRITE (ISTDE, *) 'Error when opening the tmp files'
         STOP
      END IF
      WRITE (K) 'MCP', NB
      RETURN
      END SUBROUTINE SETTMPGG
