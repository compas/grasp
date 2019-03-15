!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION SHIELD (J)
!                                                                      *
!   This  routine estimates the screening (or shielding) for orbital   *
!   J according to the relativistic hydrogenic energy formula.         *
!                                                                      *
!   Written by Yu Zou, at Vanderbilt University: 02 Mar 2000           *
!                                                                      *
!   Formula:                                                           *
!    epsilon = E(n,k)/Z^2_{eff}                                        *
!    beta = alpha Z_{eff}                                              *
!    epsilon = 1/beta^{2} [(1+beta^2/nu^2)^{-1/2} - 1]                 *
!    nu = n + (k^2 - beta^2)^{1/2} - |k|                               *
!   where E(n,k) is the orbital energy for the specific principal      *
!   quantum n and k quantum number, alpha is the fine structure        *
!   constant.                                                          *
!   Note that SHIELD is set to -C if Zeff is out of range. This will   *
!   ensure that the self-energy for such orbital is zero.              *
!   Reference: M J Seaton, Rep. Prog. Phys., Vol.46, P167,1983         *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: A, AA, BB, CC, T, BETA
!-----------------------------------------------
!
!
! a=[(1+\alpha^2 E)^{-2} - 1]^{-1/2}
! note E(J) = -E(n,k)
      A = 1 - E(J)/(C*C)
      A = 1/(A*A) - 1
      IF (A < 0.0) THEN
         SHIELD = -C
         RETURN
      ENDIF
      A = 1/SQRT(A)
!
      AA = 1 + A*A
      BB = A*DBLE(NP(J)-IABS(NAK(J)))
      CC = DBLE((NP(J)-IABS(NAK(J)))**2-NAK(J)**2)
      T = BB*BB - AA*CC
      IF (T < 0.0) THEN
         SHIELD = -C
         RETURN
      ENDIF
      BETA = BB + SQRT(T)
      BETA = BETA/AA
      SHIELD = Z - BETA*C
!     print *, j, '  Initial Screen Factor = ',shield
      RETURN
      END FUNCTION SHIELD
