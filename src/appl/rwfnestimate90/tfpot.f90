

!***********************************************************************
!                                                                      *
      SUBROUTINE TFPOT
!                                                                      *
!   Calculation of the universal Thomas-Fermi potential.               *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!                                         Last revision: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEBUG_C
      USE DEF_C, ONLY: Z, NELEC
      USE GRID_C
      USE NPOT_C, ONLY: ZZ
      USE ORB_C
      USE POTE_C, ONLY: YP
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE draw_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NB3, NROWS, II, II1, II2, II3
      REAL(DOUBLE) :: THIRD, WA, WB, WC, WD, WE, WF
!-----------------------------------------------
!
!
      THIRD = 1.0D00/3.0D00
!
      WA = Z - DBLE(NELEC - 1)
      WB = MAX(Z - WA,0.0D00)
      WB = WB**THIRD/0.8853D00
      DO I = 1, N
!
!   Rational function approximation to the universal Thomas-Fermi
!   function
!
         WC = SQRT(WB*R(I))
         WD = WC*(0.60112D0*WC + 1.81061D0) + 1.0D00
         WE = WC*(WC*(WC*(WC*(0.04793D0*WC + 0.21465D0) + 0.77112D0) + &
            1.39515D0) + 1.81061D0) + 1.0D00
         WF = WD/WE
         YP(I) = (ZZ(I)-WA)*WF*WF + WA
      END DO
!
!   Debug printout
!
      IF (LDBPR(26)) THEN
         WRITE (99, 300)
         NB3 = N/3
         IF (3*NB3 == N) THEN
            NROWS = NB3
         ELSE
            NROWS = NB3 + 1
         ENDIF
         DO II = 1, NROWS
            II1 = II
            II2 = II1 + NROWS
            II3 = II2 + NROWS
            IF (II3 <= N) THEN
               WRITE (99, 301) R(II1), YP(II1), R(II2), YP(II2), R(II3), YP(II3&
                  )
            ELSE IF (II2 <= N) THEN
               WRITE (99, 301) R(II1), YP(II1), R(II2), YP(II2)
            ELSE
               WRITE (99, 301) R(II1), YP(II1)
            ENDIF
         END DO
         CALL DRAW (YP, 1.0D00, YP, 1.0D00, N)
      ENDIF
!
      RETURN
!
  300 FORMAT(/,/,/,' Thomas-Fermi potential'/,/,3(&
         ' --------- r --------- ------ -r*V(r) ------'))
  301 FORMAT(1P,6(1X,1D21.14))
      RETURN
!
      END SUBROUTINE TFPOT
