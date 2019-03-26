!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ3(J, K, L, M, N, ITIK, SI)
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
      USE CONS_C,          ONLY: ZERO, ONE, TWO, THREE, FOUR, SEVEN
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J, K, L, M, N
      INTEGER, INTENT(IN) :: ITIK
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
         IF (IXJTIK(J,K,L,M,N,6) == 0) RETURN
      ENDIF
      I1 = (J + K + L)/2
      AS = DBLE(I1)
      A = DBLE(L)
      B = DBLE(J)
      C = DBLE(K)
      AKA = ONE
      IF (MOD(I1,2) /= 0) AKA = -AKA
      IF (J - N == 6) THEN
! -3
         IF (K - M == 6) THEN
!                        -3  -3
            SI = AKA*DSQRT((AS + ONE)*AS*(AS - ONE)*(AS - TWO)*(AS - THREE)*(AS&
                - FOUR)/((B + ONE)*B*(B - ONE)*(B - TWO)*(B - THREE)*(B - FOUR)&
               *(B - THREE - TWO)))
            SI = SI*DSQRT((AS - A - THREE - TWO)*(AS - A - FOUR)*(AS - A - &
               THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C + ONE)*C*(C - &
               ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO)))
         ELSE IF (M - K == 6) THEN
!                        3  -3
            SI = AKA*DSQRT((AS - C)*(AS - C - ONE)*(AS - C - TWO)*(AS - C - &
               THREE)*(AS - C - FOUR)*(AS - C - THREE - TWO)/((B + ONE)*B*(B - &
               ONE)*(B - TWO)*(B - THREE)*(B - FOUR)*(B - THREE - TWO)))
            SI = SI*DSQRT((AS - B + THREE + THREE)*(AS - B + THREE + TWO)*(AS&
                - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C&
                + FOUR + THREE)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR&
               )*(C + THREE)*(C + TWO)*(C + ONE)))
         ELSE IF (K - M == 4) THEN
!                       -2  -3
            SI = AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - THREE)*(AS&
                - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE - TWO)*(B - FOUR)*(&
               B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
            SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
                - A - ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)&
               *C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 4) THEN
!                        2  -3
            SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - C - FOUR)*(AS - &
               C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE)*C))
            SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - THREE - TWO)*(B - &
               FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
         ELSE IF (K - M == 2) THEN
!                       -1   -3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B&
                + TWO)*(AS - B + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B&
                - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                + ONE)))
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + &
               THREE)))
         ELSE IF (M - K == 2) THEN
!                        1   -3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS&
                - A)*(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((&
               B - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                + ONE)))
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
                - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*&
               (C + ONE)*C*(C - ONE)))
         ELSE IF (M - K == 0) THEN
!                        0   -3
            SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A&
                - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE - TWO)*(B - FOUR)*(B&
                - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)))
            SI = SI*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + FOUR)*(C + THREE)*(C&
                + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO)))
         ENDIF
      ELSE IF (N - J == 6) THEN
!  3
         IF (K - M == 6) THEN
!                        -3  3
            SI = AKA*DSQRT((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - B - &
               THREE)*(AS - B - FOUR)*(AS - B - THREE - TWO)/((C + ONE)*C*(C - &
               ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO)))
            SI = SI*DSQRT((AS - C + THREE + THREE)*(AS - C + THREE + TWO)*(AS&
                - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B&
                + FOUR + THREE)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR&
               )*(B + THREE)*(B + TWO)*(B + ONE)))
         ELSE IF (M - K == 6) THEN
!                        3   3
            SI = AKA*DSQRT((AS - A + THREE + THREE)*(AS - A + THREE + TWO)*(AS&
                - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)/((B&
                + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS + SEVEN)*(AS + THREE + THREE)*(AS + THREE + TWO)*&
               (AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + SEVEN)*(C + THREE + &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE)))
         ELSE IF (K - M == 4) THEN
!                       -2   3
            SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - B - FOUR)&
               *(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + &
               SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - FOUR)*(C - THREE)*(C&
                - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)))
         ELSE IF (M - K == 4) THEN
!                        2   3
            SI = -AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - A + THREE&
                + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
                + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + &
               FOUR)*(B + THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
               AS + THREE)*(AS + TWO)/((C + THREE + THREE)*(C + THREE + TWO)*(C&
                + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C))
         ELSE IF (K - M == 2) THEN
!                       -1   3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(&
               AS + THREE)*(AS + TWO)*(AS - B - THREE)*(AS - B - TWO)*(AS - B&
                - ONE)*(AS - B)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + &
               TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
                - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
               TWO)*(C + THREE)))
         ELSE IF (M - K == 2) THEN
!                        1   3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C&
                + TWO)*(AS - C + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A&
                + TWO)*(AS - A + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*&
               C*(C - ONE)))
         ELSE IF (M - K == 0) THEN
!                        0   3
            SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - B&
                - TWO)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + &
               SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
               THREE)*(B + TWO)*(B + ONE)))
            SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS + &
               FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(C + THREE)*(C + TWO)*&
               (C + ONE)*C*(C - ONE)*(C - TWO)))
         ENDIF
      ELSE IF (J - N == 4) THEN
! -2
         IF (K - M == 6) THEN
!                       -3  -2
            SI = AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - THREE)*(AS&
                - TWO)*(AS - ONE)*AS*(AS + ONE)/((C - THREE - TWO)*(C - FOUR)*(&
               C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
            SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
                - A - ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
               *B*(B + ONE)*(B + TWO)))
         ELSE IF (M - K == 6) THEN
!                       3   -2
            SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - C - FOUR)&
               *(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + &
               SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + &
               THREE)*(C + TWO)*(C + ONE)))
            SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
               THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - FOUR)*(B - THREE)*(B&
                - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)))
         ELSE IF (K - M == 4) THEN
!                      -2  -2
            SI = AKA*((TWO + THREE)*(AS - C)*(AS - B) - (AS + TWO)*(AS - A - &
               FOUR))*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - FOUR)*(B&
                - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
               TWO)))
         ELSE IF (M - K == 4) THEN
!                       2  -2
            SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - C + ONE)&
               *(AS - B + THREE + TWO))*DSQRT((AS - C - THREE)*(AS - C - TWO)*(&
               AS - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - &
               ONE)*B*(B + ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
                - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C&
                + THREE)*(C + TWO)*(C + ONE)*C))
         ELSE IF (K - M == 2) THEN
!                       -1  -2
            SI = AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS + TWO)*(AS - A - THREE)&
               )*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - ONE)*AS*&
               (AS + ONE)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + &
               ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE)*(&
               C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)))
         ELSE IF (M - K == 2) THEN
!                        1  -2
            SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - C + ONE)*(AS - B + FOUR))*&
               DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - C - TWO)*(AS&
                - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
               *B*(B + ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - &
               ONE)))
         ELSE IF (M - K == 0) THEN
!                        0  -2
            SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - C + ONE)*(AS - B + &
               THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - &
               ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                + ONE)*(B + TWO)))
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)&
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - &
               TWO)))
         ENDIF
      ELSE IF (N - J == 4) THEN
!  2
         IF (K - M == 6) THEN
!                        -3  2
            SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - B - FOUR)*(AS - &
               B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + &
               THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + &
               ONE)*B))
            SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - THREE - TWO)*(C - &
               FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
         ELSE IF (M - K == 6) THEN
!                        3  2
            SI = -AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - A + THREE&
                + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
                + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + &
               FOUR)*(C + THREE)*(C + TWO)*(C + ONE)))
            SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
               AS + THREE)*(AS + TWO)/((B + THREE + THREE)*(B + THREE + TWO)*(B&
                + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B))
         ELSE IF (K - M == 4) THEN
!                       -2  2
            SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - B + ONE)&
               *(AS - C + THREE + TWO))*DSQRT((AS - B - THREE)*(AS - B - TWO)*(&
               AS - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - &
               ONE)*C*(C + ONE)*(C + TWO)))
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
                - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B&
                + THREE)*(B + TWO)*(B + ONE)*B))
         ELSE IF (M - K == 4) THEN
!                        2  2
            SI = AKA*((TWO + THREE)*(AS - B)*(AS - C) - (AS - A)*(AS + THREE + &
               THREE))*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(&
               AS - A + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*&
               (B + THREE)*(B + TWO)*(B + ONE)*B))
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*&
               (C + TWO)*(C + ONE)*C))
         ELSE IF (K - M == 2) THEN
!                       -1  2
            SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - B + ONE)*(AS - C + &
               FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - B&
                - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + THREE)*(B + THREE&
                + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B))
            SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - &
               THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)))
         ELSE IF (M - K == 2) THEN
!                        1  2
            SI = -AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE + TWO&
               ))*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - A + &
               THREE)*(AS - A + TWO)*(AS - A + ONE)/((B + THREE + THREE)*(B + &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B))
            SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + THREE + TWO&
               )*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)))
         ELSE IF (M - K == 0) THEN
!                        0  2
            SI = AKA*((AS - B - TWO)*(AS - C) - (AS - A)*(AS + FOUR))*DSQRT(TWO&
               *THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS&
                - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B&
                + THREE)*(B + TWO)*(B + ONE)*B))
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - &
               TWO)))
         ENDIF
      ELSE IF (J - N == 2) THEN
! - 1
         IF (K - M == 6) THEN
!                       -3   -1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C&
                + TWO)*(AS - C + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((C&
                - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C&
                + ONE)))
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
               A)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + &
               THREE)))
         ELSE IF (M - K == 6) THEN
!                       3    -1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(&
               AS + THREE)*(AS + TWO)*(AS - C - THREE)*(AS - C - TWO)*(AS - C&
                - ONE)*(AS - C)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE + &
               TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)))
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
                - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + &
               TWO)*(B + THREE)))
         ELSE IF (K - M == 4) THEN
!                       -2   -1
            SI = AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS + TWO)*(AS - A - THREE)&
               )*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - ONE)*AS*&
               (AS + ONE)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + &
               ONE)*(C + TWO)))
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE)*(&
               B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)))
         ELSE IF (M - K == 4) THEN
!                       2    -1
            SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - C + ONE)*(AS - B + &
               FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - C&
                - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + THREE)*(C + THREE&
                + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C))
            SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - &
               THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)))
         ELSE IF (K - M == 2) THEN
!                        -1   -1
            SI = AKA*(TWO*THREE*(AS - C)*(AS - C - ONE)*(AS - B - ONE)*(AS - B)&
                - TWO*FOUR*(AS - C)*(AS - B)*(AS + TWO)*(AS - A - TWO) + (AS + &
               THREE)*(AS + TWO)*(AS - A - THREE)*(AS - A - TWO))
            SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - THREE)*(B&
                - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)*(C - THREE)*&
               (C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)))
         ELSE IF (M - K == 2) THEN
!                         1   -1
            SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) - &
               TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - C + ONE)*(AS - B + &
               THREE) + (AS - C + ONE)*(AS - C + TWO)*(AS - B + FOUR)*(AS - B&
                + THREE))
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)&
               /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + &
               THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE)*C*(C - ONE)))
         ELSE IF (M - K == 0) THEN
!                         0   -1
            SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*&
               AS*(AS - A + ONE)*(AS - C + ONE)*(AS - B + TWO) + (AS - C + ONE)&
               *(AS - C + TWO)*(AS - B + THREE)*(AS - B + TWO))
            SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)/((B&
                - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)*&
               (C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO)&
               ))
         ENDIF
      ELSE IF (N - J == 2) THEN
! - 1
         IF (K - M == 6) THEN
!                        -3   1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS&
                - A)*(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((&
               C - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C&
                + ONE)))
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
                - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*&
               (B + ONE)*B*(B - ONE)))
         ELSE IF (M - K == 6) THEN
!                        3    1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B&
                + TWO)*(AS - B + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A&
                + TWO)*(AS - A + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C + &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)))
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
               )/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*&
               B*(B - ONE)))
         ELSE IF (K - M == 4) THEN
!                        -2   1
            SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - B + ONE)*(AS - C + FOUR))*&
               DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - B - TWO)*(AS&
                - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)&
               *C*(C + ONE)*(C + TWO)))
            SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + &
               THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - &
               ONE)))
         ELSE IF (M - K == 4) THEN
!                         2   1
            SI = -AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE + TWO&
               ))*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - A + &
               THREE)*(AS - A + TWO)*(AS - A + ONE)/((C + THREE + THREE)*(C + &
               THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C))
            SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + THREE + TWO&
               )*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)))
         ELSE IF (M - K == 2) THEN
!                         1   1
            SI = AKA*(TWO*THREE*(AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE)&
                - TWO*FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + FOUR) + (AS - A)*(&
               AS - A - ONE)*(AS + THREE + TWO)*(AS + FOUR))
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B&
               *(B - ONE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C&
                + ONE)*C*(C - ONE)))
         ELSE IF (K - M == 2) THEN
!                        -1   1
            SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) - &
               TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - B + ONE)*(AS - C + &
               THREE) + (AS - B + ONE)*(AS - B + TWO)*(AS - C + FOUR)*(AS - C&
                + THREE))
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)&
               /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + &
               THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + &
               ONE)*B*(B - ONE)))
         ELSE IF (M - K == 0) THEN
!                         0   1
            SI = -AKA*TWO*((AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C - ONE&
               ) - THREE*(AS - B - ONE)*(AS - C)*(AS - A)*(AS + THREE) + (AS - &
               A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE))
            SI = SI*DSQRT(THREE*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS + &
               TWO)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + &
               ONE)*B*(B - ONE)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C&
                - ONE)*(C - TWO)))
         ENDIF
      ELSE IF (J - N == 0) THEN
!  0
         IF (K - M == 6) THEN
!                       -3    0
            SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A&
                - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE - TWO)*(C - FOUR)*(C&
                - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)))
            SI = SI*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C + &
               THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + FOUR)*(B + THREE)*(B&
                + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO)))
         ELSE IF (M - K == 6) THEN
!                        3    0
            SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - C&
                - TWO)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + &
               SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + &
               THREE)*(C + TWO)*(C + ONE)))
            SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS + &
               FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(B + THREE)*(B + TWO)*&
               (B + ONE)*B*(B - ONE)*(B - TWO)))
         ELSE IF (K - M == 4) THEN
!                       -2    0
            SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - B + ONE)*(AS - C + &
               THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - &
               ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C&
                + ONE)*(C + TWO)))
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)&
               /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B - &
               TWO)))
         ELSE IF (M - K == 4) THEN
!                        2    0
            SI = AKA*((AS - C - TWO)*(AS - B) - (AS - A)*(AS + FOUR))*DSQRT(TWO&
               *THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS&
                - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C&
                + THREE)*(C + TWO)*(C + ONE)*C))
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)&
               /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B - &
               TWO)))
         ELSE IF (K - M == 2) THEN
!                       -1    0
            SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*&
               AS*(AS - A + ONE)*(AS - B + ONE)*(AS - C + TWO) + (AS - B + ONE)&
               *(AS - B + TWO)*(AS - C + THREE)*(AS - C + TWO))
            SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)/((C&
                - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)*&
               (B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO)&
               ))
         ELSE IF (M - K == 2) THEN
!                        1    0
            SI = -AKA*TWO*((AS - C - ONE)*(AS - C - TWO)*(AS - B)*(AS - B - ONE&
               ) - THREE*(AS - C - ONE)*(AS - B)*(AS - A)*(AS + THREE) + (AS - &
               A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE))
            SI = SI*DSQRT(THREE*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS + &
               TWO)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
               ONE)*C*(C - ONE)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B&
                - ONE)*(B - TWO)))
         ELSE IF (K - M == 0) THEN
!                        0    0
            SI = AKA*((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C&
                - ONE)*(AS - C - TWO) - THREE*THREE*(AS - B)*(AS - B - ONE)*(AS&
                - C)*(AS - C - ONE)*(AS - A)*(AS + TWO) + THREE*THREE*(AS - B)*&
               (AS - C)*(AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO) - (AS&
                - A)*(AS - A - ONE)*(AS - A - TWO)*(AS + FOUR)*(AS + THREE)*(AS&
                + TWO))
            SI = SI/DSQRT((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE&
               )*(B - TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - &
               ONE)*(C - TWO))
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SIXJ3
