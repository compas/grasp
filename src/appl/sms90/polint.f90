!************************************************************************
!
      SUBROUTINE POLINT(XA, YA, N, X, Y, DY)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(DOUBLE), INTENT(IN) :: X
      REAL(DOUBLE), INTENT(OUT) :: Y, DY
      REAL(DOUBLE), DIMENSION(N),  INTENT(IN) :: XA, YA
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NMAX = 10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS, I, M
      REAL(DOUBLE), DIMENSION(NMAX) :: C, D
      REAL(DOUBLE) :: DIF, DIFT, HO, HP, W, DEN
!-----------------------------------------------
      NS = 1
      DIF = ABS(X - XA(1))
      DO I = 1, N
         DIFT = ABS(X - XA(I))
         IF (DIFT < DIF) THEN
            NS = I
            DIF = DIFT
         ENDIF
         C(I) = YA(I)
         D(I) = YA(I)
      END DO
      Y = YA(NS)
      NS = NS - 1
      DO M = 1, N - 1
         DO I = 1, N - M
            HO = XA(I) - X
            HP = XA(I+M) - X
            W = C(I+1) - D(I)
            DEN = HO - HP
            IF (DEN == 0.0D00) THEN
               WRITE (*, '(2A)') 'PAUSE ', 'FAILURE IN POLINT'
               READ *
            ENDIF
            DEN = W/DEN
            D(I) = HP*DEN
            C(I) = HO*DEN
         END DO
         IF (2*NS < N - M) THEN
            DY = C(NS+1)
         ELSE
            DY = D(NS)
            NS = NS - 1
         ENDIF
         Y = Y + DY
      END DO
      RETURN
      END SUBROUTINE POLINT
