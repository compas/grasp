!***********************************************************************
!                                                                      *
      SUBROUTINE XPOT(J)
!                                                                      *
!   This subroutine tabulates the exchange terms (the first terms on   *
!   the  right-hand  sides of eqs (14), I P Grant, B J McKenzie, P H   *
!   Norrington, D F Mayers, and N C Pyper, Computer  Phys  Commun 21   *
!   (1980) 211) for orbital J. The exchange terms are stored  in the   *
!   common arrays XP and XQ.                                           *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:59:40   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB
      USE debug_C
      USE def_C
      USE grid_C
      USE orb_C
      USE pote_C
      USE scf_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE yzk_I
      USE draw_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INDEX, LABEL, K, IOY1, IOY2, IORB, NB2, NROWS, II, II1, II2
      REAL(DOUBLE) :: COEFF, CTB
!-----------------------------------------------
!
!   Debug printout: header
!
      IF (LDBPR(27) .OR. LDBPR(28)) WRITE (99, 300) NP(J), NH(J)
!
!   Clear for accumulation of sums
!
      XP(:N) = 0.D0
      XQ(:N) = 0.D0
!
!   Add contributions from exchange terms
!
      DO INDEX = 1, NXCOF

         ! Decode information in label
         LABEL = NXA(INDEX)
         K = MOD(LABEL,KEY)
         LABEL = LABEL/KEY
         IOY1 = MOD(LABEL,KEY)
         LABEL = LABEL/KEY
         IOY2 = MOD(LABEL,KEY)
         IORB = LABEL/KEY
         COEFF = XA(INDEX)

         ! Debug printout: composition
         IF (LDBPR(27)) WRITE (99, 301) K, COEFF, NP(IOY1), NH(IOY1), NP(IOY2)&
            , NH(IOY2), NP(IORB), NH(IORB)

         CALL YZK (K, IOY1, IOY2)
!
!   Accumulate contributions
!
         COEFF = COEFF/C
         !DO I = 1, MF(IORB)
         DO I = 1, N
            CTB = COEFF*TB(I)
            XP(I) = XP(I) + CTB*QF(I,IORB)
            XQ(I) = XQ(I) - CTB*PF(I,IORB)
         END DO
      END DO
!
!   Debug printout: potential functions
!
      IF (LDBPR(28)) THEN
         WRITE (99, 302)
         NB2 = N/2
         IF (2*NB2 == N) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2 + 1
         ENDIF
         DO II = 1, NROWS
            II1 = II
            II2 = II1 + NROWS
            IF (II2 <= N) THEN
               WRITE (99, 303) R(II1), XP(II1), XQ(II1), R(II2), XP(II2), XQ(&
                  II2)
            ELSE IF (II1 <= N) THEN
               WRITE (99, 303) R(II1), XP(II1), XQ(II1)
            ENDIF
         END DO
         CALL DRAW (XP, 1.0D00, XQ, C, N)
      ENDIF
!
      RETURN
!
  300 FORMAT(/,/,' Exchange potential contributions (coefficients will ',&
         ' be divided by C) for ',1I2,1A2,' orbital :'/,/)
  301 FORMAT(/,25X,'(',1I2,')'/,1X,1P,D21.14,'* Y    (',1I2,1A2,',',1I2,1A2,&
         ') ','* P (',1I2,1A2,')')
  302 FORMAT(/,/,31X,'(P)',19X,'(Q)',41X,'(P)',19X,'(Q)'/,2(&
         ' --------- r --------- ------ X  (r) -------',&
         ' ------ X  (r) -------'))
  303 FORMAT(1P,6(1X,1D21.14))
      RETURN
!
      END SUBROUTINE XPOT
