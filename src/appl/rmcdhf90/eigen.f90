!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION EIGEN (J)
!                                                                      *
!   This function computes an estimate of the energy of orbital J .    *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE grid_C
      USE orb_C
      USE pote_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      USE dpbdt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: PIECE1, PIECE2, PIECE3, PIECE4, PIECE5
!-----------------------------------------------
!
!
!   Initialization
!
      MTP = MF(J)
!
!   Exchange term
!
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = (PF(I,J)*XQ(I)-QF(I,J)*XP(I))*RPOR(I)
      END DO
      CALL QUAD (PIECE1)
      PIECE1 = C*PIECE1
!
!   Direct term
!
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = (PF(I,J)**2 + QF(I,J)**2)*YP(I)*RPOR(I)
      END DO
      CALL QUAD (PIECE2)
!
!   Kinetic energy terms
!
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = QF(I,J)**2*RP(I)
      END DO
      CALL QUAD (PIECE3)
      PIECE3 = 2.0D00*C*C*PIECE3
!
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = (PF(I,J)*QF(I,J))*RPOR(I)
      END DO
      CALL QUAD (PIECE4)
      PIECE4 = -2.0D00*DBLE(NAK(J))*C*PIECE4
!
      CALL DPBDT (J)
      TA(1) = 0.0D00
      DO I = 2, MTP
         TA(I) = PF(I,J)*TB(I) - QF(I,J)*TA(I)
      END DO
      CALL QUAD (PIECE5)
      PIECE5 = C*PIECE5/H
!
!   Assembly
!
      EIGEN = PIECE1 + PIECE2 + PIECE3 + PIECE4 + PIECE5
!
      RETURN
!
      END FUNCTION EIGEN
