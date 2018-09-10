!*******************************************************************
!                                                                  *
      SUBROUTINE SUWJJ(K1,K2,LL,J1,J2,SUW)
!                                                                  *
!                 (k1 k2)                                          *
!     ( j QJ ::: W      ::: j QJ )                                 *
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
      USE CONS_C,          ONLY: ZERO, EPS
      USE ribojj_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumtjj_I
      USE ixjtik_I
      USE rmeajj_I
      USE sixj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: K1, K2, J1, J2
      INTEGER,      INTENT(OUT) :: LL
      REAL(DOUBLE), INTENT(OUT) :: SUW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK1,KK2,I,IP,IG,L2Q1,L2V1,L2J1,L2Q2,L2V2,L2J2, &
                 L2QI,L2VI,L2JI
      REAL(DOUBLE) :: COEF1, COEF2, S, SI1, SI2
!-----------------------------------------------
      SUW=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      S=ZERO
      CALL RUMTJJ(J1,LL,L2Q1,L2V1,L2J1)
      KK1=K1*2
      KK2=K2*2
      IP=IMPNJJ(J1)
      IG=IMGNJJ(J1)
      CALL RUMTJJ(J2,LL,L2Q2,L2V2,L2J2)
      DO I=IP,IG
        CALL RUMTJJ(I,LL,L2QI,L2VI,L2JI)
        IF(IXJTIK(LL,LL,KK2,L2J2,L2J1,L2JI) /= 0) THEN
          IF(IXJTIK(1,1,KK1,L2Q2,L2Q1,L2QI) /= 0) THEN
            CALL RMEAJJ(LL,J1,L2Q1,L2J1,I,L2QI,L2JI,COEF1)
            IF(DABS(COEF1) > EPS) THEN
              CALL RMEAJJ(LL,I,L2QI,L2JI,J2,L2Q2,L2J2,COEF2)
              IF(DABS(COEF2) > EPS) THEN
                CALL SIXJ(LL,LL,KK2,L2J2,L2J1,L2JI,0,SI1)
                CALL SIXJ(1,1,KK1,L2Q2,L2Q1,L2QI,0,SI2)
                S=S+SI1*SI2*COEF1*COEF2
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      END DO
      SUW=S*DSQRT(DBLE((KK1+1)*(KK2+1)))
      IF(MOD(L2Q1+L2J1+L2Q2+L2J2+KK1+KK2,4) /= 0)SUW=-SUW
      RETURN
      END SUBROUTINE SUWJJ
