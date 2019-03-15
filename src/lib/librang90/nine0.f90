!*******************************************************************
!                                                                  *
      SUBROUTINE NINE0(J1,J2,J3,L1,L2,L3,K1,K2,K3,AA)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
!                                                                  *
!     |  J1/2  J2/2  J3/2 |                                        *
!     |  L1/2  L2/2  L3/2 |                                        *
!     |  K1/2  K2/2    0  |                                        *
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
      USE CONS_C,          ONLY: ZERO, ONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J1, J2, J3, L1, L2, L3, K1, K2, K3
      REAL(DOUBLE), INTENT(OUT) :: AA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IFA
      REAL(DOUBLE) :: A, B
!-----------------------------------------------
      IF (J1 == 0) THEN
         CALL SIXJ (L2, K2, J2, K3, L3, L1, 0, A)
         B = DBLE((J2 + 1)*(L1 + 1))
         IFA = K2 + J2 + L3 + L1
      ELSE IF (J2 == 0) THEN
         CALL SIXJ (L3, K3, J3, K1, L1, L2, 0, A)
         B = DBLE((J3 + 1)*(L2 + 1))
         IFA = K3 + J3 + L1 + L2
      ELSE IF (J3 == 0) THEN
         CALL SIXJ (L1, K1, J1, K2, L2, L3, 0, A)
         B = DBLE((J1 + 1)*(L3 + 1))
         IFA = K1 + J1 + L2 + L3
      ELSE IF (L1 == 0) THEN
         CALL SIXJ (K2, J2, L2, J3, K3, K1, 0, A)
         B = DBLE((L2 + 1)*(K1 + 1))
         IFA = J2 + L2 + K3 + K1
      ELSE IF (L2 == 0) THEN
         CALL SIXJ (K3, J3, L3, J1, K1, K2, 0, A)
         B = DBLE((L3 + 1)*(K2 + 1))
         IFA = J3 + L3 + K1 + K2
      ELSE IF (L3 == 0) THEN
         CALL SIXJ (K1, J1, L1, J2, K2, K3, 0, A)
         B = DBLE((L1 + 1)*(K3 + 1))
         IFA = J1 + L1 + K2 + K3
      ELSE IF (K1 == 0) THEN
         CALL SIXJ (J2, J3, J1, L3, L2, K2, 0, A)
         B = DBLE((J1 + 1)*(K2 + 1))
         IFA = J3 + J1 + L2 + K2
      ELSE IF (K2 == 0) THEN
         CALL SIXJ (J3, J1, J2, L1, L3, K3, 0, A)
         B = DBLE((J2 + 1)*(K3 + 1))
         IFA = J1 + J2 + L3 + K3
      ELSE IF (K3 == 0) THEN
         CALL SIXJ (J1, J2, J3, L2, L1, K1, 0, A)
         B = DBLE((J3 + 1)*(K1 + 1))
         IFA = J2 + J3 + L1 + K1
      ELSE
         A = ZERO
         B = ONE
      ENDIF
      AA = A/DSQRT(B)
      IF (MOD(IFA,4) /= 0) AA = -AA
      RETURN
      END SUBROUTINE NINE0
