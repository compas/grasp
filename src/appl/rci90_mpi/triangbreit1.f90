!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TRIANGBREIT1 (IA,IB,IC,ID,K)
!                                                                      *
!     Evaluate the triangular relation for the Breit integral of       *
!     See paper by I. Grant and B McKenzie, JPB 13, (1980) 2671-2681   *
!     formula 5                                                        *
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
      INTEGER, INTENT(IN) :: IA,IB,IC,ID,K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LA,LB,LC,LD,JA,JB,JC,JD,L
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

      TRIANGBREIT1 = .FALSE.

      DO L = K-1,K+1
         IF ((MOD(K+LA+LB,2).EQ.1).AND.(ABS(JA-JB).LE.2*L).AND.     &
            (ABS(JA+JB).GE.2*L)) THEN
            IF ((MOD(K+LC+LD,2).EQ.1).AND.(ABS(JC-JD).LE.2*L).AND.  &
               (ABS(JC+JD).GE.2*L)) THEN
               TRIANGBREIT1 = .TRUE.
               RETURN
            END IF
         END IF
      END DO

      RETURN
      END FUNCTION TRIANGBREIT1
