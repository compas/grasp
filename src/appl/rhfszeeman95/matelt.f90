!***********************************************************************
!                                                                      *
      SUBROUTINE MATELT(I1, K, I2, APART, GJPART, DGJPART) 
!                                                                      *
!   This  routine computes  the angular part  of the  reduced matrix   *
!   elements of the magnetic and electric multipole operators as       *
!   well as for th operators associated with the g_j factor            *
!                                                                      *
!   Calls to: [LIB92]: CLRX                                            *
!                                                                      *
!   Written by Per Jonsson                Last revision: 22 Oct 1999   *
!                                                                      *
!   Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03 1/3/07*  
!   Modified by Charlotte Froese Fischer                               *
!                     Gediminas Gaigalas  11/01/17                     *
!   Modified by Wenxian Li F77 to F90 12/28/18                         *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW
      USE orb_C,           ONLY: nak
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE clrx_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I1 
      INTEGER  :: K 
      INTEGER, INTENT(IN) :: I2 
      REAL(DOUBLE), INTENT(OUT) :: APART 
      REAL(DOUBLE), INTENT(OUT) :: GJPART 
      REAL(DOUBLE), INTENT(OUT) :: DGJPART 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: KAP1, KAP2, L1, L2 
      REAL(DOUBLE) :: FASE, OVLFAC 
!-----------------------------------------------
!
!
!   Set KAP1 and KAP2
!
      KAP1 = NAK(I1) 
!
      IF (MOD(K,2) == 1) THEN 
         KAP2 = -NAK(I2) 
      ELSE 
         KAP2 = NAK(I2) 
      ENDIF 
!
!   Determine the l quantum numbers
!
      IF (KAP1 > 0) THEN 
         L1 = KAP1 
      ELSE 
         L1 = (-KAP1) - 1 
      ENDIF 
!
      IF (KAP2 > 0) THEN 
         L2 = KAP2 
      ELSE 
         L2 = (-KAP2) - 1 
      ENDIF 
!
      IF (MOD(L1 + K + L2,2) == 0) THEN 
!
!   Parity selection rule satisfied
!
!   Determine the phase factor
!
         IF (MOD(KAP1 + 1,2) == 0) THEN 
            FASE = 1.0D00 
         ELSE 
            FASE = -1.0D00 
         ENDIF 
!
!   The other factor is \sqrt (2 j_2 + 1); since j = | \kappa | - 1/2,
!   we have 2 j_2 + 1 = 2 | \kappa |; the factor \sqrt (2 j_2 + 1)
!   has been accounted for in MCT
!
         OVLFAC = FASE*SQRT(DBLE(2*ABS(KAP2))) 
!
         IF (MOD(K,2) == 1) THEN 
!
!   These are for the magnetic multipole moments and the two operators of the g_j factor
!
            APART = (KAP1 + NAK(I2))*CLRX(KAP1,K,KAP2)*OVLFAC 
!
            GJPART = APART 
            DGJPART = -(KAP1 + NAK(I2)-1.0D00)*CLRX(KAP1,K,KAP2)*OVLFAC 
         ELSE 
!
!   These are for the electric multipole moments
!
            APART = CLRX(KAP1,K,KAP2)*OVLFAC 
!
         ENDIF 
!
      ELSE 
!
         APART = 0.0D00 
         GJPART = 0.0D00 
         DGJPART = 0.0D00 
!
      ENDIF 
!
      RETURN  
      END SUBROUTINE MATELT 
