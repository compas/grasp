!*******************************************************************
!                                                                  *
      SUBROUTINE COULOM(J1,J2,J3,J4,L1,L2,L3,L4,K,AA)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 01  ------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
!                                                                  *
!                          k   k+1  (k) (k)                        *
!     (n l j T  n l j T ::r  / r  ( C   C )::n l j T  n l j T )    *
!       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
!                                                                  *
!     SUBROUTINE CALLED:  CRE                                      *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C,          ONLY:  EPS, ZERO
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      USE cre_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: J1,J2,J3,J4,L1,L2,L3,L4,K
      REAL(DOUBLE), INTENT(OUT)  :: AA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IFAZ
!-----------------------------------------------
      AA=ZERO
      IF(ITTK(L1,L3,K) == 0)RETURN
      IF(ITTK(L2,L4,K) == 0)RETURN
      I=(2*K+1)/2
      AA=CRE (J1,I,J3)
      IF(DABS(AA) < EPS)RETURN
      AA=AA*CRE (J2,I,J4)
      IF(DABS(AA) < EPS)RETURN
      IFAZ=L3-2*K-L1+L4-L2
      IF((IFAZ/4)*4 /= IFAZ)AA=-AA
      RETURN
      END SUBROUTINE COULOM
