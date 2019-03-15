!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION BRINTF (ITYPE, IA, IB, IC, ID, K)
!                                                                      *
!   Computes integrals for the transverse photon interaction.          *
!                                                                      *
!   Call(s) to: [RCI92]: BESSEL, BRRA.                                 *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def,   ONLY: KEYORB
      USE stor_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE bessel_I
      USE brra_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ITYPE
      INTEGER  :: IA
      INTEGER  :: IB
      INTEGER  :: IC
      INTEGER  :: ID
      INTEGER  :: K
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: KEY2 = KEY*KEY
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOD, ICOD1, ICOD2, ICOD3, ICOD4
!-----------------------------------------------
!
      GO TO (1,4,3,6,2,5) ITYPE
!
!   Type 1 and 5 integrals require j(k), n(k) Bessel fuctions
!   Type 5 integrals only require w = wab Bessel functions
!
    1 CONTINUE
      IF (IA/=IC .OR. IB/=ID .OR. IA/=ID .OR. IC==IB) CALL BESSEL (IC, ID, 1, 2&
         , K)
    2 CONTINUE
      CALL BESSEL (IA, IB, 1, 1, K)
      GO TO 6
!
!   Type 3 integrals require j(k), n(k) Bessel functions for either
!   w = wab or w = cd whichever is non-zero.
!
    3 CONTINUE
      IF (IA /= IB) CALL BESSEL (IA, IB, 1, 1, K)
      IF (IC /= ID) CALL BESSEL (IC, ID, 1, 2, K)
      GO TO 6
!
!   Type 2 and 6 integrals require j(k), n(k) and j(k+2), n(k+2)
!   Bessel fuctions
!   Type 6 integrals only require w = wab Bessel functions.
!
    4 CONTINUE
      IF (IA/=IC .OR. IB/=ID .OR. IA/=ID .OR. IC/=IB) THEN
!
         ICOD = MAX(IC,ID) + KEY*MIN(IC,ID)
         ICOD1 = ICOD + KEY2*(K - 1)
         ICOD2 = ICOD + KEY2*(K + 1)
         ICOD = MAX(IA,IB) + KEY*MIN(IA,IB)
         ICOD3 = ICOD + KEY2*(K - 1)
         ICOD4 = ICOD + KEY2*(K + 1)
         IF (ICOD1==KEEP(1,2) .AND. ICOD2==KEEP(2,2) .AND. ICOD3==KEEP(1,1)&
             .AND. ICOD4==KEEP(2,1)) GO TO 6
         IF (ICOD1==KEEP(1,1) .AND. ICOD2==KEEP(2,1) .AND. ICOD3==KEEP(1,2)&
             .AND. ICOD4==KEEP(2,2)) GO TO 6
         CALL BESSEL (IC, ID, 1, 2, K - 1)
         CALL BESSEL (IC, ID, 2, 2, K + 1)
      ENDIF
!
    5 CONTINUE
      CALL BESSEL (IA, IB, 1, 1, K - 1)
      CALL BESSEL (IA, IB, 2, 1, K + 1)
!
!   Compute the integral
!
    6 CONTINUE
      BRINTF = BRRA(ITYPE,IA,IB,IC,ID,K)
!
      RETURN
      END FUNCTION BRINTF
