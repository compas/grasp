!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION SKINT (RAC, IA, IC, RBD, IB, ID, K, IW)
!                                                                      *
!   This routine evaluates transverse interaction integrals:           *
!                                                                      *
!                          (k)                                         *
!                         S   (a,c;b,d;w)                              *
!                                                                      *
!   where w = wac if IW = 1, and w= wbd if IW = 2.                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 06 Nov 1992   *
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
      USE grid_C, ONLY: r, rp, rpor
      USE orb_C
      USE tatb_C, ONLY: mtp, ta, tb
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
      INTEGER  :: IA, IB, IC, ID
      INTEGER, INTENT(IN) :: K, IW
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: RAC
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: RBD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MXRAC, MXRBD, I
      REAL(DOUBLE), DIMENSION(NNNP) :: TKEEP
      REAL(DOUBLE) :: EPSI, W, WK, VALU
!-----------------------------------------------
!
!
      EPSI = 1.0D-10
!
      W = WIJ(IW)
      WK = DBLE(K + K + 1)
      MXRAC = MIN(MF(IA),MF(IC))
      MXRBD = MIN(MF(IB),MF(ID))
!
!           (k-1)
!  Compute Z     (rho  ; s)
!                    ac
!
      TA(:MXRAC) = RAC(:MXRAC)
      MTP = MXRAC
      CALL ZKF (K - 1, IA, IC)
!
      IF (DABS(W) < EPSI) THEN
!
!   W = 0 case
!
         TKEEP(:MTP) = TB(:MTP)
         CALL ZKF (K + 1, IA, IC)
         MTP = MIN(MTP,MXRBD)
         TA(1) = 0.0D00
         TA(2:MTP) = RBD(2:MTP)*RPOR(2:MTP)*(TB(2:MTP)-TKEEP(2:MTP))
         CALL QUAD (VALU)
         SKINT = WK*VALU*0.5D00
!
      ELSE
!
!   Finite w: see I P Grant and B J McKenzie, J Phys B: At Mol Phys,
!   13 (1980) 2671-2681
!
         TKEEP(:MTP) = TB(:MTP)
         TA(:MTP) = -TA(:MTP)*BESSJ(1,IW,:MTP)
         CALL ZKF (K - 1, IA, IC)
!
         MTP = MIN(MTP,MXRBD)
         TA(1) = 0.0D00
         TA(2:MTP) = ((1.0D00 + BESSN(2,IW,2:MTP))*TB(2:MTP)-TKEEP(2:MTP)*BESSN&
            (2,IW,2:MTP))*RBD(2:MTP)/R(2:MTP)**2*RPOR(2:MTP)
         CALL QUAD (VALU)
         SKINT = (WK/W)**2*VALU
         TA(:MXRBD) = RBD(:MXRBD)*(1.0D00 + BESSJ(2,IW,:MXRBD))
         MTP = MXRBD
         CALL ZKF (K + 1, IB, ID)
         MTP = MIN(MTP,MXRAC)
         TA(1) = 0.0D00
         TA(2:MTP) = RAC(2:MTP)*(1.0D00 + BESSN(1,IW,2:MTP))*TB(2:MTP)*R(2:MTP)&
            *RP(2:MTP)
         CALL QUAD (VALU)
         SKINT = SKINT - VALU*W*W/DBLE((2*K + 3)*(2*K - 1))
!
      ENDIF
!
      IF (LDBPR(11)) WRITE (99, 300) K, NP(IA), NH(IA), NP(IC), NH(IC), NP(IB)&
         , NH(IB), NP(ID), NH(ID), IW, VALU
!
      RETURN
!
  300 FORMAT('  (',1I2,')'/,'S     (',1I2,1A2,',',1I2,1A2,'|',1I2,1A2,',',1I2,&
         1A2,';',1I2,') = ',1P,D19.12)
      RETURN
!
      END FUNCTION SKINT
