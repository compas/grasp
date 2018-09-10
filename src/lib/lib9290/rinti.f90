!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINTI (J, K, MODE) 
!                                                                      *
!   The value of this  function is the one-electron integral I (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
!   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:23   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE DEBUG_C
      USE DEF_C,           ONLY: C 
      USE GRID_C 
      USE NPOT_C,          ONLY: ZZ 
      USE ORB_C 
      USE TATB_C,          ONLY: TA, TB, MTP
      USE WAVE_C,          ONLY: MF,PF,QF 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dpbdt_I 
      USE quad_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
      INTEGER  :: K 
      INTEGER, INTENT(IN) :: MODE 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: PIECE1, PIECE2, PIECE3, PIECE4 
!-----------------------------------------------
!
!   Stop if orbitals J and K have different kappa values
!
      IF (NAK(J) /= NAK(K)) THEN 
         WRITE (*, 300) NP(J), NH(J), NP(K), NH(K) 
         STOP  
      ENDIF 
!
      MTP = MAX(MF(J),MF(K)) 
!
!   Kinetic energy contribution
!
!   Piece involving derivatives
!
      CALL DPBDT (K) 
      TA(1) = 0.D0 
      DO I = 2, MTP 
         TA(I) = QF(I,J)*TA(I) - PF(I,J)*TB(I) 
      END DO 
      CALL QUAD (PIECE1) 
      PIECE1 = C*PIECE1/H 
!
!   Pieces not involving derivatives
!
      TA(1) = 0.D0 
      DO I = 2, MTP 
         TA(I) = RP(I)*QF(I,J)*QF(I,K) 
      END DO 
      CALL QUAD (PIECE2) 
      PIECE2 = -2.D0*C*C*PIECE2 
!
      TA(1) = 0.D0 
      DO I = 2, MTP 
         TA(I) = RPOR(I)*(PF(I,J)*QF(I,K) + QF(I,J)*PF(I,K)) 
      END DO 
      CALL QUAD (PIECE3) 
      PIECE3 = PIECE3*C*DBLE(NAK(K)) 
!
!   Contribution from nuclear potential only if MODE is 0
!
      IF (MODE == 0) THEN 
!
         TA(1) = 0.D0 
         DO I = 2, MTP 
            TA(I) = RPOR(I)*ZZ(I)*(PF(I,J)*PF(I,K) + QF(I,J)*QF(I,K)) 
         END DO 
         CALL QUAD (PIECE4) 
         PIECE4 = -PIECE4 
!
      ELSE 
         PIECE4 = 0.D0 
      ENDIF 
!
      RINTI = PIECE1 + PIECE2 + PIECE3 + PIECE4 
!
!   Debug printout
!
      IF (MODE==0 .AND. LDBPR(4)) WRITE (99, 301) NP(J), NH(J), NP(K), NH(K), &
         RINTI 
      IF (MODE/=0 .AND. LDBPR(5)) WRITE (99, 302) NP(J), NH(J), NP(K), NH(K), &
         RINTI 
!
      RETURN  
!
  300 FORMAT('RINTI: Attempt to calculate I(',1I2,1A2,',',1I2,1A2,')') 
  301 FORMAT(/,' I (',1I2,1A2,',',1I2,1A2,') = ',1P,D19.12,/) 
  302 FORMAT(/,' K (',1I2,1A2,',',1I2,1A2,') = ',1P,D19.12,/) 
      RETURN  
!
      END FUNCTION RINTI 
