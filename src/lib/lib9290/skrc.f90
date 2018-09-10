!***********************************************************************
!                                                                      *
      SUBROUTINE SKRC(IS, KAPS, KS, KD1, KD2, KE1, KE2) 
!                                                                      *
!   Determines the range of the tensor rank k for Coulomb integral.    *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:43   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: KD1 
      INTEGER, INTENT(OUT) :: KD2 
      INTEGER, INTENT(OUT) :: KE1 
      INTEGER, INTENT(OUT) :: KE2 
      INTEGER, DIMENSION(4), INTENT(IN) :: IS
      INTEGER, DIMENSION(4), INTENT(IN) :: KAPS
      INTEGER, DIMENSION(4), INTENT(IN) :: KS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISD1, ISD2, KD1A, KD1B, KD2A, KD2B, ISE1, ISE2, KE1A, &
         KE1B, KE2A, KE2B 
!-----------------------------------------------
!
!
      KD2 = 0 
      KE2 = 0 
!
!   Direct terms --- KD1 = minimum k , KD2 = number of terms
!
      ISD1 = 1 
      IF (KAPS(1)*KAPS(3) < 0) ISD1 = -1 
      ISD2 = 1 
      IF (KAPS(2)*KAPS(4) < 0) ISD2 = -1 
      KD1A = ABS(KS(1)-KS(3)) 
      IF (ISD1 < 0) KD1A = KD1A + 2 
      KD1B = ABS(KS(2)-KS(4)) 
      IF (ISD2 < 0) KD1B = KD1B + 2 
      IF (MOD((KD1A - KD1B)/2,2) == 0) THEN 
         KD2A = KS(1) + KS(3) - 2 
         IF (ISD1 > 0) KD2A = KD2A - 2 
         KD2B = KS(2) + KS(4) - 2 
         IF (ISD2 > 0) KD2B = KD2B - 2 
         KD1 = MAX(KD1A,KD1B)/2 
         KD2 = MIN(KD2A,KD2B)/2 
         KD2 = (KD2 - KD1)/2 + 1 
      ENDIF 
!
!   Exchange terms --- KE1 = minimum k , KE2 = number of terms
!
      IF (IS(1)==IS(2) .OR. IS(3)==IS(4)) RETURN  
      ISE1 = 1 
      IF (KAPS(1)*KAPS(4) < 0) ISE1 = -1 
      ISE2 = 1 
      IF (KAPS(2)*KAPS(3) < 0) ISE2 = -1 
      KE1A = ABS(KS(1)-KS(4)) 
      IF (ISE1 < 0) KE1A = KE1A + 2 
      KE1B = ABS(KS(2)-KS(3)) 
      IF (ISE2 < 0) KE1B = KE1B + 2 
      IF (MOD((KE1A - KE1B)/2,2) /= 0) RETURN  
      KE2A = KS(1) + KS(4) - 2 
      IF (ISE1 > 0) KE2A = KE2A - 2 
      KE2B = KS(2) + KS(3) - 2 
      IF (ISE2 > 0) KE2B = KE2B - 2 
      KE1 = MAX(KE1A,KE1B)/2 
      KE2 = MIN(KE2A,KE2B)/2 
      KE2 = (KE2 - KE1)/2 + 1 
!
      RETURN  
      END SUBROUTINE SKRC 
