!***********************************************************************
!                                                                      *
      SUBROUTINE SET_CSF_number
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
      USE rang_Int_C,     only: NUM_in_BLK
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: IOS, I_num, I_Num_BLK
      CHARACTER :: RECORD*1, S_closed*1, S_orbitals*2
!-----------------------------------------------
!
      NUM_in_BLK(1:20) = 0
      I_num = 0
      I_Num_BLK = 1
      DO
         READ (21, '(A2)', IOSTAT=IOS) S_orbitals
         IF (IOS == 0) THEN
            IF (S_orbitals(1:2) == ' *') THEN
               NUM_in_BLK(I_Num_BLK) = I_num
               I_Num_BLK = I_Num_BLK + 1
               I_num = 0
               READ (21, '(A2)') S_orbitals
            ENDIF
            I_num = I_num + 1
!     Read the J_sub and v quantum numbers
            READ (21, '(A2)') S_orbitals
!     Read the X, J, and (sign of) P quantum numbers
            READ (21, '(A2)') S_orbitals
            CYCLE 
         END IF
         EXIT
      END DO
      NUM_in_BLK(I_Num_BLK) = I_num
      I_Num_BLK = I_Num_BLK + 1
!---------------------------------------------------
      I_Num_BLK = 1
      REWIND(21)
      READ (21, '(1A1)') RECORD
!     Closed orbitals
      READ (21, '(A1)') S_closed
      READ (21, '(1A1)') RECORD
!     Peel orbitals
      READ (21, '(A2)') S_orbitals
      READ (21, '(1A1)') RECORD
      RETURN  
      END SUBROUTINE SET_CSF_number
