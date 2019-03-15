!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION VINTI (J, K)
!                                                                      *
!   The value of this  function is the one-electron integral V (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (3.23) in  F A  Parpia, M. Tong and C F Fischer,   *
!   to appear.                                                         *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by M Tong and F A Parpia,     Last revision: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE debug_C
      USE grid_C
      USE orb_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dpbdt_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      INTEGER  :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, KPJ, KPK, IFACT1, IFACT2
      REAL(DOUBLE) :: PIECE1, FACT1, FACT2, PIECE2
!-----------------------------------------------
!
      MTP = MAX(MF(J),MF(K))
!
!   Piece involving derivatives
!
      CALL DPBDT (K)
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = PF(I,J)*TA(I) + QF(I,J)*TB(I)
      END DO
      CALL QUAD (PIECE1)
      PIECE1 = PIECE1/H
!
!   Pieces not involving derivatives
!
      KPJ = NAK(J)
      KPK = NAK(K)
      IFACT1 = KPJ*(KPJ + 1) - KPK*(KPK + 1)
      FACT1 = 0.5D00*DBLE(IFACT1)
      IFACT2 = (-KPJ*((-KPJ) + 1)) + KPK*((-KPK) + 1)
      FACT2 = 0.5D00*DBLE(IFACT2)
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = RPOR(I)*(FACT1*PF(I,J)*PF(I,K) + FACT2*QF(I,J)*QF(I,K))
      END DO
      CALL QUAD (PIECE2)
!
      VINTI = PIECE1 - PIECE2
!
!   Debug printout
!
      IF (LDBPR(6)) WRITE (99, 300) NP(J), NH(J), NP(K), NH(K), VINTI
!
      RETURN
!
  300 FORMAT(/,'VINTI: V (',1I2,1A2,',',1I2,1A2,') = ',1P,D19.12)
      RETURN
!
      END FUNCTION VINTI
