!
!     ------------------------------------------------------------------
!     B N D I N V
!     ------------------------------------------------------------------
!
      SUBROUTINE BNDINV(A, EL, N, DETERM, EPSIL, ITEST, NSIZE)
!
!       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
!       FROM "DLYTAP".
!
!*      DOUBLE PRECISION E,F
!*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: ITEST
      INTEGER , INTENT(IN) :: NSIZE
      REAL(DOUBLE) , INTENT(OUT) :: DETERM
      REAL(DOUBLE) , INTENT(IN) :: EPSIL
      REAL(DOUBLE) , INTENT(INOUT) :: A(NSIZE,NSIZE)
      REAL(DOUBLE) , INTENT(INOUT) :: EL(NSIZE,NSIZE)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISL2, K000FX, INDSNL, I, J, N1, M, K, J1, I1, KS
      REAL(DOUBLE) :: D, C, S, DETERP, F, E, EPSILP, RAT
!-----------------------------------------------
      IF (N < 2) GO TO 140
      ISL2 = 0
      K000FX = 2
      IF (ISL2 == 0) INDSNL = 2
      IF (ISL2 == 1) INDSNL = 1
!       CALL SLITET(2,INDSNL)
!       CALL OVERFL(K000FX)
!       CALL DVCHK(K000FX)
!
!       SET EL = IDENTITY MATRIX
      DO I = 1, N
         EL(I,:N) = 0.0D0
         EL(I,I) = 1.0D0
      END DO
!
!       TRIANGULARIZE A, FORM EL
!
      N1 = N - 1
      M = 2
      DO J = 1, N1
         DO I = M, N
            IF (A(I,J) == 0.0D0) CYCLE
            D = DSQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
            C = A(J,J)/D
            S = A(I,J)/D
            DO K = J, N
               D = C*A(J,K) + S*A(I,K)
               A(I,K) = C*A(I,K) - S*A(J,K)
               A(J,K) = D
            END DO
            DO K = 1, N
               D = C*EL(J,K) + S*EL(I,K)
               EL(I,K) = C*EL(I,K) - S*EL(J,K)
               EL(J,K) = D
            END DO
         END DO
         M = M + 1
      END DO
!       CALL OVERFL(K000FX)
!       GO TO (140,51),K000FX
!
!       CALCULATE THE DETERMINANT
      DETERP = A(1,1)
      DO I = 2, N
         DETERP = DETERP*A(I,I)
      END DO
      DETERM = DETERP
!       CALL OVERFL(K000FX)
!       GO TO (140,520,520),K000FX
!
!       IS MATRIX SINGULAR
      F = A(1,1)
      E = A(1,1)
      DO I = 2, N
         IF (DABS(F) < DABS(A(I,I))) F = A(I,I)
         IF (DABS(E) <= DABS(A(I,I))) CYCLE
         E = A(I,I)
      END DO
      EPSILP = EPSIL
      IF (EPSILP <= 0) EPSILP = 1.0E-8
      RAT = E/F
      IF (ABS(RAT) < EPSILP) GO TO 130
!
!       INVERT TRIANGULAR MATRIX
      J = N
      DO J1 = 1, N
!       CALL SLITE(2)
         I = J
         ISL2 = 1
         DO I1 = 1, J
!       CALL SLITET(2,K000FX)
            IF (ISL2 == 0) K000FX = 2
            IF (ISL2 == 1) THEN
               K000FX = 1
               ISL2 = 0
            ENDIF
            SELECT CASE (K000FX)
            CASE DEFAULT
               A(I,J) = 1.0D0/A(I,I)
            CASE (2)
               KS = I + 1
               D = 0.0D0
               D = SUM(A(I,KS:J)*A(KS:J,J))
               A(I,J) = -D/A(I,I)
            END SELECT
 1003       CONTINUE
            I = I - 1
         END DO
         J = J - 1
      END DO
!       CALL OVERFL(K000FX)
!       GO TO (140,103,103),K000FX

!103    CALL DVCHK(K000FX)
!       GO TO (140,105),K000FX
!
!       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
      M = 1
      DO I = 1, N
         DO J = 1, N
            D = 0.0D0
            D = SUM(A(I,M:N)*EL(M:N,J))
            EL(I,J) = D
         END DO
         M = M + 1
      END DO
!       CALL OVERFL(K000FX)
!       GO TO (140,123,123),K000FX
!
!       RECOPY EL TO A
      A(:N,:N) = EL(:N,:N)
      ITEST = 0
!126    IF(INDSNL.EQ.1)CALL SLITE(2)
  126 CONTINUE
      IF (INDSNL == 1) ISL2 = 1
      RETURN
!
  130 CONTINUE
      ITEST = 1
      GO TO 126
  140 CONTINUE
      ITEST = -1
      GO TO 126
      END SUBROUTINE BNDINV
