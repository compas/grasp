!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ5(J,K,L,M,N,ITIK,SI)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  1/2 |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE, TWO
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER                   :: J, K, L, M, N
      INTEGER, INTENT(IN)       :: ITIK
      REAL(DOUBLE), INTENT(OUT) :: SI
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: I1
      REAL(DOUBLE) :: AS, A, B, C, AKA
!-----------------------------------------------
      SI = ZERO
      IF (ITIK /= 0) THEN
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(J,K,L,M,N,1) == 0) RETURN
      ENDIF
      I1 = (J + K + L)/2
      AS = DBLE(I1)
      A = DBLE(L)
      B = DBLE(K)
      C = DBLE(J)
      AKA = ONE
      IF (MOD(I1,2) /= 0) AKA = -AKA
      IF (K < M) THEN
         IF (J < N) THEN
!              M > K,  J < N.
            SI = -AKA*DSQRT((AS + TWO)*(AS - A + ONE)/((B + ONE)*(B + TWO)*(C&
                + ONE)*(C + TWO)))
         ELSE IF (J > N) THEN
!              M > K,  J > N.
            SI = AKA*DSQRT((AS - C + ONE)*(AS - B)/((B + ONE)*(B + TWO)*C*(C + &
               ONE)))
         ENDIF
      ELSE IF (K > M) THEN
         IF (J < N) THEN
!             M < K,  J < N.
            SI = AKA*DSQRT((AS - C)*(AS - B + ONE)/(B*(B + ONE)*(C + ONE)*(C + &
               TWO)))
         ELSE IF (J > N) THEN
!             M < K,  J > N.
            SI = AKA*DSQRT((AS + ONE)*(AS - A)/(B*(B + ONE)*C*(C + ONE)))
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SIXJ5
