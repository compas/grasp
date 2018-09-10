!***********************************************************************
!                                                                      *
      SUBROUTINE LAGCON(J, NPROCS) 
!                                                                      *
!   This  routine  includes  the Lagrange multiplier contribution in   *
!   the 'exchange' term.                                               *
!
!   Parameter nprocs is added so that it works for mpi and serial
!   programs.
!                                                                      *
!                                           Last update: 08 Dec 1992   *
!   Modified by Xinghong He                 Last update: 17 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: KEYORB
      USE def_C
      USE grid_C
      USE lagr_C
      USE pote_C
      USE scf_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J 
      INTEGER , INTENT(IN) :: NPROCS 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IECCK, L1, L2, M, MFM, I 
      REAL(DOUBLE) :: EPS, WB, WA, WARI 
!-----------------------------------------------
!
!
      IF (NEC == 0) RETURN  
!
!ww      EPS = ACCY*0.1D 00
      EPS = ACCY*0.01D00 
!
!   Add contributions from off-diagonal parameters to exchange
!
      WB = 1.0D00/(UCF(J)*C)/NPROCS 
      DO K = 1, NEC 
!
!   Decode index
!
         IECCK = IECC(K) 
         L1 = IECCK/KEY 
         L2 = IECCK - KEY*L1 
!
         IF (J == L1) THEN 
            M = L2 
         ELSE IF (J == L2) THEN 
            M = L1 
         ELSE 
            CYCLE  
         ENDIF 
!
         WA = ECV(K)*WB 
         IF (ABS(WA) < EPS) CYCLE  
!
!   ADD CONTRIBUTIONS TO EXCHANGE TERMS
!
         MFM = MF(M) 
         DO I = 1, MFM 
            WARI = WA*R(I) 
            XP(I) = XP(I) + WARI*QF(I,M) 
            XQ(I) = XQ(I) - WARI*PF(I,M) 
         END DO 
!
      END DO 
!
      RETURN  
      END SUBROUTINE LAGCON 
