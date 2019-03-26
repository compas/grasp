!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RKINT (RAC, IA, IC, RBD, IB, ID, K, IW)
!                                                                      *
!   This routine evaluates the transverse interaction integrals.  If   *
!   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
!   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
!   w = wbd  if IW = 2.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 15 Oct 1992   *
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
      USE bess_C
      USE debug_C
      USE grid_C, ONLY: rpor
      USE orb_C
      USE tatb_C,  ONLY: mtp, ta, tb
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE zkf_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA, IC
      INTEGER, INTENT(IN) :: IB, ID
      INTEGER  :: K
      INTEGER, INTENT(IN) :: IW
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: RAC
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: RBD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MXRBD, MXRAC, I
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------
!
      MXRBD = MIN (MF(IB),MF(ID))
      MXRAC = MIN (MF(IA),MF(IC))
!
      IF (IW == 0) THEN
!
!   IW = 0
!
         TA(:MXRAC) = RAC(:MXRAC)
         MTP = MXRAC
         CALL ZKF (K, IA, IC)
         MTP = MIN(MTP,MXRBD)
         TA(1) = 0.0D00
         TA(2:MTP) = RBD(2:MTP)*TB(2:MTP)*RPOR(2:MTP)
!
      ELSE
!
!   IW = 1,2
!
         TA(:MXRAC) = RAC(:MXRAC)*(1.0D00 + BESSJ(1,IW,:MXRAC))
         MTP = MXRAC
         CALL ZKF (K, IA, IC)
         MTP = MIN(MTP,MXRBD)
         TA(1) = 0.0D00
         TA(2:MTP) = RBD(2:MTP)*(1.0D00 + BESSN(1,IW,2:MTP))*TB(2:MTP)*RPOR(2:&
            MTP)
!
      ENDIF
!
      CALL QUAD (RESULT)
      RKINT = RESULT
!
!   Debug printout if option set
!
      IF (LDBPR(11)) WRITE (99, 300) K, NP(IA), NH(IA), NP(IC), NH(IC), NP(IB)&
         , NH(IB), NP(ID), NH(ID), IW, RESULT
!
      RETURN
!
  300 FORMAT('_ (',1I2,')'/,'R     (',1I2,1A2,',',1I2,1A2,'|',1I2,1A2,',',1I2,&
         1A2,';',1I2,') = ',1P,D19.12)
      RETURN
!
      END FUNCTION RKINT
