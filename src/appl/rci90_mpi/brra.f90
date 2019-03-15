!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION BRRA (ITYPE, IA, IC, IB, ID, K)
!                                                                      *
!   This routine evaluates the transverse interaction integrals:       *
!                                                                      *
!      ITYPE = 1: General R (k; a c | b d )                            *
!            = 2: General S (k; a c | b d )                            *
!            = 3: R (k; a , b d )                                      *
!            = 4: F (k; a , b )                                        *
!            = 5: G (k; a , b ) type integral                          *
!            = 6: H (k; a , b ) type integral                          *
!                                                                      *
!   Call(s) to: [RCI92]: RKINT, SKINT.                                 *
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
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE cons_C
      USE grid_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rkint_I
      USE skint_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ITYPE
      INTEGER  :: IA
      INTEGER  :: IC
      INTEGER  :: IB
      INTEGER  :: ID
      INTEGER  :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MXRAC, I, MXRBD
      REAL(DOUBLE), DIMENSION(NNNP) :: RAC, RBD
!-----------------------------------------------
!
      MXRAC = MIN(MF(IA),MF(IC))
      DO I = 1, MXRAC
         RAC(I) = PF(I,IA)*QF(I,IC)
      END DO
!
      MXRBD = MIN(MF(IB),MF(ID))
      DO I = 1, MXRBD
         RBD(I) = PF(I,IB)*QF(I,ID)
      END DO
!
      GO TO (21,22,23,24,25,26) ITYPE
!
!   ITYPE = 1
!
   21 CONTINUE
      IF (IA==IB .AND. IC==ID) GO TO 9
      IF (IA==ID .AND. IC==IB) GO TO 10
      BRRA = (RKINT(RAC,IA,IC,RBD,IB,ID,K,1) + RKINT(RAC,IA,IC,RBD,IB,ID,K,2)&
          + RKINT(RBD,IB,ID,RAC,IA,IC,K,1) + RKINT(RBD,IB,ID,RAC,IA,IC,K,2))*&
         HALF
      RETURN
!
!   ITYPE = 2
!
   22 CONTINUE
      IF (IA==IB .AND. IC==ID) GO TO 26
      IF (IA==ID .AND. IC==IB) GO TO 26
      BRRA = (SKINT(RAC,IA,IC,RBD,IB,ID,K,1) + SKINT(RAC,IA,IC,RBD,IB,ID,K,2))*&
         HALF
      RETURN
!
!   ITYPE = 3
!
   23 CONTINUE
      IF (IA == IC) THEN
         DO I = 1, MXRBD
            RBD(I) = RBD(I) + PF(I,ID)*QF(I,IB)
         END DO
         BRRA = (RKINT(RAC,IA,IC,RBD,IB,ID,K,0) + RKINT(RBD,IB,ID,RAC,IA,IC,K,0&
            ) + RKINT(RAC,IA,IC,RBD,IB,ID,K,2) + RKINT(RBD,IB,ID,RAC,IA,IC,K,2)&
            )*HALF
         RETURN
      ENDIF
      DO I = 1, MXRAC
         RAC(I) = RAC(I) + PF(I,IC)*QF(I,IA)
      END DO
      BRRA = (RKINT(RAC,IA,IC,RBD,IB,ID,K,1) + RKINT(RBD,IB,ID,RAC,IA,IC,K,1)&
          + RKINT(RAC,IA,IC,RBD,IB,ID,K,0) + RKINT(RBD,IB,ID,RAC,IA,IC,K,0))*&
         HALF
      RETURN
!
!   ITYPE = 4
!
   24 CONTINUE
      BRRA = RKINT(RAC,IA,IC,RBD,IB,ID,K,0) + RKINT(RBD,IB,ID,RAC,IA,IC,K,0)
      RETURN
!
!   ITYPE = 5
!
   25 CONTINUE
      IF (IA==ID .AND. IC==IB) GO TO 10
    9 CONTINUE
      BRRA = TWO*RKINT(RAC,IA,IC,RBD,IB,ID,K,1)
      RETURN
   10 CONTINUE
      BRRA = RKINT(RAC,IA,IC,RBD,IB,ID,K,1) + RKINT(RBD,IB,ID,RAC,IA,IC,K,1)
      RETURN
!
!   ITYPE = 6
!
   26 CONTINUE
      BRRA = SKINT(RAC,IA,IC,RBD,IB,ID,K,1)
      RETURN
!
      END FUNCTION BRRA
