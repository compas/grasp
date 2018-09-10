!*******************************************************************
!                                                                  *
      SUBROUTINE RMEAJJ11(J1,J2,LL,S)
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, HALF, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumtjj_I
      USE c0t5s_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: LL, J1, J2
      REAL(DOUBLE), INTENT(OUT) :: S
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: LQ, LV, LQS, LVS, J, J1S, N, INN
      REAL(DOUBLE) :: A1, A4, Q, QQ, QS, QM, QMS
!-----------------------------------------------
      S=ZERO
      CALL RUMTJJ(J1,LL,LQ,LV,J)
      CALL RUMTJJ(J2,LL,LQS,LVS,J1S)
      QQ=HALF
      N=2
      Q=HALF*DBLE(LQ)
      QS=HALF*DBLE(LQS)
      QM=-HALF*DBLE(HALF*(LL+1)-N)
      QMS=-HALF*DBLE(HALF*(LL+1)-(N-1))
      CALL C0T5S(QS,QMS,QQ,Q,QM,A4)
      IF(DABS(A4) < EPS) RETURN
      A1=DBLE(N*(LQ+1)*(J+1))
      S=DSQRT(A1)/A4
      INN=-N-1
      IF((INN/2)*2 /= INN)S=-S
      RETURN
      END SUBROUTINE RMEAJJ11
