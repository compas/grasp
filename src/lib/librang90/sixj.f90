!*******************************************************************
!                                                                  *
      SUBROUTINE SIXJ(I,J,K,L,M,N,ITIK,SI)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
!                                                                  *
!     | I/2  J/2  K/2 |                                            *
!     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
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
      USE CONS_C,          ONLY: ZERO, ONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      USE sixj5_I
      USE sixj1_I
      USE sixj35_I
      USE sixj2_I
!      USE gracah1_I
      USE dracah_I
      USE sixj3_I
      USE sixj4_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER             :: I, J, K, L, M, N
      INTEGER, INTENT(IN) :: ITIK
      REAL(DOUBLE)        :: SI
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IFA
      REAL(DOUBLE), DIMENSION(0:4,0:4,0:4,0:4,0:4,0:4) :: RACA
      REAL(DOUBLE) :: UNDEF, A
      LOGICAL :: SAVE
!-----------------------------------------------
      DATA RACA/ 15625*1.D-20/
      DATA UNDEF/ 1.D-20/
      SI = ZERO
      IF (ITIK /= 0) THEN
!
!     CHESKED TRIANGULAR CONDITIONS
!
         IF (IXJTIK(I,J,K,L,M,N) == 0) RETURN
      ENDIF
      SAVE = .FALSE.
      IF (MAX0(I,J,K,L,M,N) <= 4) THEN
         SI = RACA(I,J,K,L,M,N)
         IF (SI == UNDEF) THEN
            SAVE = .TRUE.
         ELSE
            RETURN
         ENDIF
      ENDIF
!
!     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
!
      IF (I*J*K*L*M*N == 0) THEN
         IF (I == 0) THEN
            A = DBLE((M + 1)*(K + 1))
            IFA = L + M + K
         ELSE IF (J == 0) THEN
            A = DBLE((L + 1)*(K + 1))
            IFA = I + M + N
         ELSE IF (K == 0) THEN
            A = DBLE((I + 1)*(L + 1))
            IFA = I + M + N
         ELSE IF (L == 0) THEN
            A = DBLE((J + 1)*(K + 1))
            IFA = I + J + K
         ELSE IF (M == 0) THEN
            A = DBLE((I + 1)*(K + 1))
            IFA = I + J + K
         ELSE
            A = DBLE((I + 1)*(J + 1))
            IFA = I + J + K
         ENDIF
         SI = ONE/DSQRT(A)
         IF (MOD(IFA,4) /= 0) SI = -SI
!
!     THE CASE 1/2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 1) THEN
         IF (I == 1) THEN
            CALL SIXJ5 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 1) THEN
            CALL SIXJ5 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 1) THEN
            CALL SIXJ5 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 1) THEN
            CALL SIXJ5 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 1) THEN
            CALL SIXJ5 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ5 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 1
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 2) THEN
         IF (I == 2) THEN
            CALL SIXJ1 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 2) THEN
            CALL SIXJ1 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 2) THEN
            CALL SIXJ1 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 2) THEN
            CALL SIXJ1 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 2) THEN
            CALL SIXJ1 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ1 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 3/2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 3) THEN
         IF (I == 3) THEN
            CALL SIXJ35 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 3) THEN
            CALL SIXJ35 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 3) THEN
            CALL SIXJ35 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 3) THEN
            CALL SIXJ35 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 3) THEN
            CALL SIXJ35 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ35 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 4) THEN
         IF (I == 4) THEN
            CALL SIXJ2 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 4) THEN
            CALL SIXJ2 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 4) THEN
            CALL SIXJ2 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 4) THEN
            CALL SIXJ2 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 4) THEN
            CALL SIXJ2 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ2 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 5/2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 5) THEN
         CALL DRACAH (I, J, M, L, K, N, SI)
         IF (MOD(I + J + M + L,4) /= 0) SI = -SI
!
!     CASES 3
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 6) THEN
         IF (I == 6) THEN
            CALL SIXJ3 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 6) THEN
            CALL SIXJ3 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 6) THEN
            CALL SIXJ3 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 6) THEN
            CALL SIXJ3 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 6) THEN
            CALL SIXJ3 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ3 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 7/2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 7) THEN
         CALL DRACAH (I, J, M, L, K, N, SI)
         IF (MOD(I + J + M + L,4) /= 0) SI = -SI
!
!     CASES 4
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 8) THEN
         IF (I == 8) THEN
            CALL SIXJ4 (M, K, L, J, N, 0, SI)
         ELSE IF (J == 8) THEN
            CALL SIXJ4 (I, N, M, L, K, 0, SI)
         ELSE IF (K == 8) THEN
            CALL SIXJ4 (I, M, N, L, J, 0, SI)
         ELSE IF (L == 8) THEN
            CALL SIXJ4 (J, K, I, M, N, 0, SI)
         ELSE IF (M == 8) THEN
            CALL SIXJ4 (I, K, J, L, N, 0, SI)
         ELSE
            CALL SIXJ4 (I, J, K, L, M, 0, SI)
         ENDIF
!
!     THE CASE 9/2
!
      ELSE IF (MIN0(I,J,K,L,M,N) == 9) THEN
         CALL DRACAH (I, J, M, L, K, N, SI)
         IF (MOD(I + J + M + L,4) /= 0) SI = -SI
!
!     CALCULATED OTHER CASES
!
      ELSE
       CALL DRACAH(I,J,M,L,K,N,SI)
!         CALL GRACAH1 (I, J, M, L, K, N, SI)
!        CALL GRACAH(I,J,M,L,K,N,SI)
         IF (MOD(I + J + M + L,4) /= 0) SI = -SI
      ENDIF
      IF (SAVE) RACA(I,J,K,L,M,N) = SI
      RETURN
      END SUBROUTINE SIXJ
