!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TRIANGBREIT2 (IA,IB,IC,ID,L)
!                                                                      *
!     Evaluate the triangle relations for Breit integrals of type 2    *
!     See Grant and McKenzier JPB 13 (1980), 2675-2681 formula 5       *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE orb_C, ONLY: NKL, NKJ
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: IA,IB,IC,ID,L
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LA,LB,LC,LD,JA,JB,JC,JD
!-----------------------------------------------
!
      LA = NKL(IA)
      LB = NKL(IB)
      LC = NKL(IC)
      LD = NKL(ID)
      JA = NKJ(IA)
      JB = NKJ(IB)
      JC = NKJ(IC)
      JD = NKJ(ID)

      TRIANGBREIT2 = .FALSE.

      IF ((MOD(L+LA+LB,2).EQ.0).AND.(ABS(JA-JB).LE.2*L).AND.           &
         (ABS(JA+JB).GE.2*L)) THEN
         IF ((MOD(L+LC+LD,2).EQ.0).AND.(ABS(JC-JD).LE.2*L).AND.        &
            (ABS(JC+JD).GE.2*L)) THEN
            TRIANGBREIT2 = .TRUE.
            RETURN
         END IF
      END IF

      RETURN
      END FUNCTION TRIANGBREIT2
