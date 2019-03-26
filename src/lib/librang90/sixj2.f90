!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ2(J,K,L,M,N,ITIK,SI)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | J/2  K/2  L/2 |                                            *
!     | M/2  N/2   2  |             [B.M.X. 75]                    *
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
      USE CONS_C,          ONLY: ZERO, HALF, ONE, TWO, THREE, FOUR
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
      REAL(DOUBLE) :: AS, A, B, C, AKA, X1, X2, X3
!-----------------------------------------------
      SI = ZERO
      IF (ITIK /= 0) THEN
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(J,K,L,M,N,4) == 0) RETURN
      ENDIF
      I1 = (J + K + L)/2
      AS = DBLE(I1)
      A = DBLE(L)
      B = DBLE(J)
      C = DBLE(K)
      AKA = ONE
      IF (MOD(I1,2) /= 0) AKA = -AKA
      IF (J - N == 4) THEN
! -2
         IF (K - M == 4) THEN
!  I                      -2  -2
            SI = AKA*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE)*(B&
                - TWO)*(B - ONE)*B*(B + ONE)))
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 4) THEN
!  V   P(12)               2  -2
            SI = AKA*DSQRT((AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS&
                - C)/((C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)*(C + TWO + &
               THREE)))
            SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)*(AS - &
               B + FOUR)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ELSE IF (K - M == 2) THEN
!  II   P(12)              -1   -2
            SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((C - TWO)*(C - ONE)*C*&
               (C + ONE)*(C + TWO)))
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - C)*(AS&
                - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ELSE IF (M - K == 2) THEN
!  IV   P(12)               1   -2
            SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - C - TWO)*(AS - C - ONE&
               )*(AS - C)/(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)))
            SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)/((B - &
               THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ELSE IF (K - M == 0) THEN
!  III  P(12)               0   -2
            SI = AKA*DSQRT(TWO*THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C&
                - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)))
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)&
               /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ENDIF
!  2
      ELSE IF (N - J == 4) THEN
         IF (K - M == 4) THEN
!  V                      -2   2
            SI = AKA*DSQRT((AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS&
                - B)/((B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)*(B + TWO + &
               THREE)))
            SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)*(AS - &
               C + FOUR)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 4) THEN
!  1                       2   2
            SI = AKA*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS&
                - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*&
               (B + ONE)))
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))&
               )
         ELSE IF (K - M == 2) THEN
!  3                      -1   2
            SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - B - TWO)*(AS - B - &
               ONE)*(AS - B)/((B + TWO + THREE)*(B + FOUR)*(B + THREE)*(B + TWO&
               )*(B + ONE)))
            SI = SI*TWO*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((&
               C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 2) THEN
!  2                       1   2
            SI = -AKA*DSQRT((AS - B)*(AS - C + ONE)*(AS - A + THREE)*(AS - A + &
               TWO)*(AS - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B&
                + TWO)*(B + ONE)))
            SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(&
               C + THREE)*(C + TWO)*(C + ONE)*C))
         ELSE IF (K - M == 0) THEN
!  5                        0   2
            SI = AKA*DSQRT(THREE*TWO*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS&
                - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*&
               (B + ONE)))
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)))
         ENDIF
      ELSE IF (J - N == 2) THEN
! -1
         IF (K - M == 4) THEN
!  II   P(12)              -2  -1
            SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((B - TWO)*(B - ONE)*B*&
               (B + ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - B)*(AS&
                - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 4) THEN
!  3   P(12)                2  -1
            SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - C - TWO)*(AS - C - &
               ONE)*(AS - C)/((C + TWO + THREE)*(C + FOUR)*(C + THREE)*(C + TWO&
               )*(C + ONE)))
            SI = SI*TWO*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((&
               B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)))
         ELSE IF (K - M == 2) THEN
!  VI                       -1  -1
            SI = AKA*((A + B)*(A - B + TWO) - (C - TWO)*(C - B + TWO))/DSQRT((B&
                - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))
            SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C - TWO)*(C&
                - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 2) THEN
!  VIII P(12)              1  -1
            SI = AKA*((A + C + FOUR)*(A - C - TWO) - (B - TWO)*(B + C + FOUR))/&
               DSQRT(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)&
               /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)))
         ELSE IF (K - M == 0) THEN
!  VII  P(12)              0  -1
            SI = AKA*HALF*((A + C + TWO)*(A - C) - B*B + FOUR)/DSQRT((C - ONE)*&
               C*(C + ONE)*(C + TWO)*(C + THREE))
            SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)&
               /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)))
         ENDIF
      ELSE IF (N - J == 2) THEN
!  1
         IF (K - M == 4) THEN
!  IV                     -2   1
            SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - B - TWO)*(AS - B - ONE&
               )*(AS - B)/(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)))
            SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)/((C - &
               THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 4) THEN
!  2 P(12)                 2   1
            SI = -AKA*DSQRT((AS - C)*(AS - B + ONE)*(AS - A + THREE)*(AS - A + &
               TWO)*(AS - A + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C&
                + TWO)*(C + ONE)))
            SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(&
               B + THREE)*(B + TWO)*(B + ONE)*B))
         ELSE IF (K - M == 2) THEN
!  VIII                   -1   1
            SI = AKA*((A + B + FOUR)*(A - B - TWO) - (C - TWO)*(B + C + FOUR))/&
               DSQRT(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)&
               /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 2) THEN
!  4                       1   1
            SI = AKA*(THREE*(AS - B)*(AS - C) - (AS - A)*(AS + FOUR))/DSQRT((B&
                + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C))
         ELSE IF (K - M == 0) THEN
!  6 P(12)                 0   1
            SI = -AKA*((AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE))/DSQRT((&
               B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)
            SI = SI*DSQRT(THREE*TWO*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS&
                + TWO)/((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)))
         ENDIF
      ELSE IF (N - J == 0) THEN
! 0
         IF (K - M == 4) THEN
!  III                     -2   0
            SI = AKA*DSQRT(THREE*TWO*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B&
                - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)))
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)&
               /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 4) THEN
!  5                        2   0
            SI = AKA*DSQRT(THREE*TWO*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS&
                - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*&
               (C + ONE)))
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)))
         ELSE IF (K - M == 2) THEN
!  VII                     -1   0
            SI = AKA*HALF*((A + B + TWO)*(A - B) - C*C + FOUR)/DSQRT((B - ONE)*&
               B*(B + ONE)*(B + TWO)*(B + THREE))
            SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)&
               /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 2) THEN
!  6                        1   0
            SI = -AKA*((AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE))/DSQRT((&
               C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)
            SI = SI*DSQRT(THREE*TWO*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS&
                + TWO)/((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)))
         ELSE IF (K - M == 0) THEN
!  IX                       0   0
            X1 = (AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE)
            X2 = FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + TWO)
            X3 = (AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO)
            SI = AKA*(X1 - X2 + X3)/DSQRT((B - ONE)*B*(B + ONE)*(B + TWO)*(B + &
               THREE)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SIXJ2
