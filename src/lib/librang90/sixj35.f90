!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ35(J,K,L,M,N,ITIK,SI)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2  3/2 |             [B.M.X. 75]                    *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE, TWO, THREE, FOUR
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
      INTEGER :: I1
      REAL(DOUBLE) :: AS, A, B, C, AKA
!-----------------------------------------------
      SI = ZERO
      IF (ITIK /= 0) THEN
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(J,K,L,M,N,3) == 0) RETURN
      ENDIF
      I1 = (J + K + L)/2
      AS = DBLE(I1)
      A = DBLE(L)
      B = DBLE(J)
      C = DBLE(K)
      AKA = ONE
      IF (MOD(I1,2) /= 0) AKA = -AKA
      IF (J - N == 3) THEN
! -3
         IF (K - M == 3) THEN
!  I                      -3/2  -3/2
            SI = AKA*DSQRT((AS - ONE)*AS*(AS + ONE)*(AS - A - TWO)*(AS - A - &
               ONE)*(AS - A)/((B - TWO)*(B - ONE)*B*(B + ONE)*(C - TWO)*(C - &
               ONE)*C*(C + ONE)))
         ELSE IF (M - K == 3) THEN
!  IV  P(12)              3/2   -3/2
            SI = AKA*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + ONE&
               )*(AS - B + TWO)*(AS - B + THREE)/((C + 1)*(C + TWO)*(C + THREE)&
               *(C + FOUR)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ELSE IF (K - M == 1) THEN
!  II  P(12)             -1/2   -3/2
            SI = AKA*DSQRT(THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - C)&
               *(AS - B + ONE)/((C - ONE)*C*(C + ONE)*(C + TWO)*(B - TWO)*(B - &
               ONE)*B*(B + ONE)))
         ELSE IF (M - K == 1) THEN
!  III P(12)              1/2   -3/2
            SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C - ONE)*(AS - C)*(&
               AS - B + ONE)*(AS - B + TWO)/(C*(C + ONE)*(C + TWO)*(C + THREE)*&
               (B - TWO)*(B - ONE)*B*(B + ONE)))
         ENDIF
      ELSE IF (N - J == 3) THEN
!  3
         IF (K - M == 3) THEN
!  IV                     -3/2   3/2
            SI = AKA*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C + ONE&
               )*(AS - C + TWO)*(AS - C + THREE)/((B + ONE)*(B + TWO)*(B + &
               THREE)*(B + FOUR)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 3) THEN
!  2       pataisyta               3/2   3/2
            SI = -AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS + FOUR)*(AS - A + ONE)*&
               (AS - A + TWO)*(AS - A + THREE)/((B + ONE)*(B + TWO)*(B + THREE)&
               *(B + FOUR)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)))
         ELSE IF (K - M == 1) THEN
!  1   P(12)   pataisytas          -1/2    3/2
            SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - C + ONE)*(AS&
                - C + TWO)*(AS - B - ONE)*(AS - B)/((C - ONE)*C*(C + ONE)*(C + &
               TWO)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)))
         ELSE IF (M - K == 1) THEN
!  3  P(12)     taisyta           1/2    3/2
            SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A&
                + TWO)*(AS - B)*(AS - C + ONE)/(C*(C + ONE)*(C + TWO)*(C + &
               THREE)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)))
         ENDIF
! -1
      ELSE IF (J - N == 1) THEN
         IF (K - M == 3) THEN
!  II                   -3/2   -1/2
            SI = AKA*DSQRT((THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - B&
               )*(AS - C + ONE))/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - TWO)*(C&
                - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 3) THEN
!  1                     3/2   -1/2
            SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - B + ONE)*(AS&
                - B + TWO)*(AS - C - ONE)*(AS - C)/((B - ONE)*B*(B + ONE)*(B + &
               TWO)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)))
         ELSE IF (K - M == 1) THEN
!  V                    -1/2   -1/2
            SI = AKA*(TWO*(AS - B)*(AS - C) - (AS + TWO)*(AS - A - ONE))*DSQRT(&
               (AS + ONE)*(AS - A)/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - ONE)*C&
               *(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 1) THEN
!  VI P(12)              1/2   -1/2
            SI = AKA*((AS - B + TWO)*(AS - C + ONE) - TWO*(AS - A + ONE)*(AS + &
               ONE))*DSQRT((AS - C)*(AS - B + ONE)/(C*(C + ONE)*(C + TWO)*(C + &
               THREE)*(B - ONE)*B*(B + ONE)*(B + TWO)))
         ENDIF
! 1
      ELSE IF (N - J == 1) THEN
         IF (K - M == 3) THEN
!  III                  -3/2    1/2
            SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B - ONE)*(AS - B)*(&
               AS - C + ONE)*(AS - C + TWO)/(B*(B + ONE)*(B + TWO)*(B + THREE)*&
               (C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 3) THEN
!  3              pataisyta       3/2    1/2
            SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A&
                + TWO)*(AS - B + ONE)*(AS - C)/(B*(B + ONE)*(B + TWO)*(B + &
               THREE)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)))
         ELSE IF (K - M == 1) THEN
!  VI                   -1/2    1/2
            SI = AKA*((AS - C + TWO)*(AS - B + ONE) - TWO*(AS - A + ONE)*(AS + &
               ONE))*DSQRT((AS - B)*(AS - C + ONE)/(B*(B + ONE)*(B + TWO)*(B + &
               THREE)*(C - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 1) THEN
!  4      pataisyta               1/2    1/2
            SI = -AKA*(TWO*(AS - B)*(AS - C) - (AS + THREE)*(AS - A))*DSQRT((AS&
                + TWO)*(AS - A + ONE)/(B*(B + ONE)*(B + TWO)*(B + THREE)*C*(C&
                + ONE)*(C + TWO)*(C + THREE)))
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SIXJ35
