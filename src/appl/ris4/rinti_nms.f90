!***********************************************************************
!                                                                      *
       SUBROUTINE RINTI_NMS (J,K,VALNMSK123,VALNMSK1)
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
!   Modified by C. Naz\'e  Oct. 2011                                   *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNN1
      USE def_C
      USE grid_C
      USE wave_C
      USE orb_C
      USE tatb_C
      USE debug_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dpbdt_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(OUT) :: VALNMSK123, VALNMSK1
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(IN) :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNN1) :: PA, PB, QA, QB
      REAL(DOUBLE) :: A1, A2, PIECE1, PIECE2, PIECE3, PIECE4
      INTEGER :: L, KAP1, KAP2, I, L1, L2, J1, J2, L2_TILDA
!-----------------------------------------------
!
!
!   Set KAP1 and KAP2
!
!C*      RINTI_NMS = 0.0D 00
      VALNMSK1   = 0.0D00
      VALNMSK123 = 0.0D00
      KAP1 = NAK(J)
      KAP2 = NAK(K)
!
!   Stop if orbitals J and K have different kappa values
!
      IF (KAP1 .NE. KAP2) THEN
         WRITE (*,300) NP(J),NH(J),NP(K),NH(K)
         STOP
      ENDIF
!
!   Determine the l quantum numbers
!
      IF (KAP1 .GT. 0) THEN
         L1 =  KAP1
      ELSE
         L1 = -KAP1-1
      ENDIF
!
      IF (KAP2 .GT. 0) THEN
         L2 =  KAP2
      ELSE
         L2 = -KAP2-1
      ENDIF
      IF (L1 .NE. L2) RETURN
!
!   Determine the j quantum numbers and l_2 tilda
!
      J1 = 2*IABS(NAK(J))-1
      J2 = 2*IABS(NAK(K))-1
      L2_TILDA = J2 - L2
!
      MTP = MAX (MF(J),MF(K))
!
!   Kinetic energy contribution
!
!   Piece involving derivatives
!
      CALL DPBDT (J)
      PA(1) = 0.0D00
      QA(1) = 0.0D00
      DO 11 I = 2,MTP
        PA(I) = TA(I)
        QA(I) = TB(I)
   11 CONTINUE
      CALL DPBDT (K)
      PB(1) = 0.0D00
      QB(1) = 0.0D00
      DO 12 I = 2,MTP
        PB(I) = TA(I)
        QB(I) = TB(I)
   12 CONTINUE
      TA(1) = 0.0D00
      DO 1 I = 2,MTP
         TA(I) = (PA(I)*PB(I)+QA(I)*QB(I))/RP(I)
    1 CONTINUE
      CALL QUAD (PIECE1)
      PIECE1 = PIECE1/(H*H)
!
      A1 = DBLE(L2*(L2+1))
      A2 = DBLE(L2_TILDA*(L2_TILDA+1))
      TA(1) = 0.0D00
      DO 2 I = 2,MTP
        TA(I) = RP(I)                                                  &
              *(A1*PF(I,J)*PF(I,K)+A2*QF(I,J)*QF(I,K))                 &
              /(R(I)*R(I))
    2 CONTINUE
      CALL QUAD (PIECE2)
!
!   Pieces not involving derivatives
!
      TA(1) = 0.D0
      DO 3 I = 2,MTP
         TA(I) = (QF(I,J)*PB(I)+QF(I,K)*PA(I))/R(I)
    3 CONTINUE
      CALL QUAD (PIECE3)
      PIECE3 = -2.0*Z*PIECE3/(C*H)
!
      IF (KAP2 .NE. 1) THEN
        TA(1) = 0.0D00
        DO 4 I = 2,MTP
           TA(I) = RP(I)                                               &
                 * (QF(I,J)*PF(I,K)+PF(I,J)*QF(I,K))                   &
                 /(R(I)*R(I))
    4   CONTINUE
        CALL QUAD (PIECE4)
        PIECE4 = -PIECE4*Z*(DBLE (KAP2)-1.0)/C
      ELSE
        PIECE4 = 0.0D00
      END IF
      VALNMSK1 = 0.5*(PIECE1+PIECE2)
      VALNMSK123 = 0.5*(PIECE1+PIECE2+PIECE3+PIECE4)
!
!   Debug printout
!
!C*      LDBPR(5)=.TRUE.
      IF (LDBPR(5))                                                    &
         WRITE (99,302) NP(J),NH(J),NP(K),NH(K),VALNMSK123
!
!C*      RETURN
!
  300 FORMAT ('RINTI_NMS: Attempt to calculate I(',                    &
               1I2,1A2,',',1I2,1A2,')')
  302 FORMAT (/' K (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
!
      RETURN
      END SUBROUTINE RINTI_NMS
