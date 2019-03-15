!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION CRE (KAP1, K, KAP2)
!-----------------------------------------------
!                                                                      *
!   Computes the relativistic reduced matrix element                   *
!                                                                      *
!                         (j1 || C(K) || j2),                          *
!                                                                      *
!   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
!   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
!   conditions are tested by the routine CLRX.                         *
!                                                                      *
!   Call(s) to: [LIB92] CLRX.                                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:10   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE clrx_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KAP1
      INTEGER  :: K
      INTEGER  :: KAP2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K1
      REAL(DOUBLE) :: DK1K2
!-----------------------------------------------
!
      K1 = ABS(KAP1)
      DK1K2 = DBLE(4*K1*IABS(KAP2))
      CRE = SQRT(DK1K2)*CLRX(KAP1,K,KAP2)
      IF (MOD(K1,2) == 1) CRE = -CRE
!
      RETURN
      END FUNCTION CRE
