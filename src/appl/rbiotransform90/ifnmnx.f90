!
!     ------------------------------------------------------------------
!     I F N M N X
!     ------------------------------------------------------------------
!
      INTEGER FUNCTION IFNMNX (IVEC, NEL, IMXMN) 
!
! Smallest or largest value of integer array
!
! IMXMN = 1 => Largest value
! IMXMN = 2 => Smallest
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NEL 
      INTEGER , INTENT(IN) :: IMXMN 
      INTEGER , INTENT(IN) :: IVEC(NEL) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IVAL, IEL, NTEST 
!-----------------------------------------------
!
      IVAL = IVEC(1) 
      IF (IMXMN == 1) THEN 
         IVAL = MAX0(MAXVAL(IVEC(2:NEL)),IVAL) 
      ELSE IF (IMXMN == 2) THEN 
         IVAL = MIN0(MINVAL(IVEC(2:NEL)),IVAL) 
      ELSE 
         WRITE (6, *) ' Stop in IFNMNX ' 
         WRITE (6, *) ' Improper calue of IMXMN ', IMXMN 
         STOP 'IFNMNX' 
      ENDIF 
!
      IFNMNX = IVAL 
!
      NTEST = 0 
      IF (NTEST /= 0) WRITE (6, *) ' Value returned from IFNMNX', IFNMNX 
!
      RETURN  
      END FUNCTION IFNMNX 
