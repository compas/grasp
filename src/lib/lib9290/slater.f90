!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION SLATER (IA, IB, IC, ID, K) 
!                                                                      *
!   The value of this  function is the Slater integral                 *
!                                                                      *
!                               k                                      *
!                              R (abcd)                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD, YZK.                                    *
!                                                                      *
!                                         Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:44   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE DEBUG_C 
      USE GRID_C 
      USE ORB_C 
      USE TATB_C,          ONLY: TA, TB, MTP
      USE WAVE_C,          ONLY: MF, QF, PF 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE yzk_I 
      USE quad_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: IA 
      INTEGER :: IB 
      INTEGER :: IC 
      INTEGER :: ID 
      INTEGER :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: RESULT 
!-----------------------------------------------
!
!
      CALL YZK (K, IB, ID) 
!
!   Multiply by second term, and obtain result by integration
!
      IF (K==0 .AND. IB==ID) THEN 
         MTP = MIN(MF(IA),MF(IC)) 
      ELSE 
         MTP = MIN(MIN(MTP,MF(IA)),MF(IC)) 
      ENDIF 
!
      TA(1) = 0.0D00 
      DO I = 2, MTP 
         TA(I) = (PF(I,IA)*PF(I,IC) + QF(I,IA)*QF(I,IC))*RPOR(I)*TB(I) 
      END DO 
!
      CALL QUAD (RESULT) 
      SLATER = RESULT 
!
!   Debug printout
!
      IF (LDBPR(10)) WRITE (99, 300) K, NP(IA), NH(IA), NP(IB), NH(IB), NP(IC)&
         , NH(IC), NP(ID), NH(ID), SLATER 
!
      RETURN  
!
  300 FORMAT(/,'  (',1I1,')',/,' R   (',1I2,1A2,',',1I2,1A2,';',1I2,1A2,',',1I2&
         ,1A2,') ','= ',1P,D19.12) 
      RETURN  
!
      END FUNCTION SLATER 
