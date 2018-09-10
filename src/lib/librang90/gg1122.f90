!*******************************************************************
!                                                                  *
      SUBROUTINE GG1122(K1,K2,QM1,QM2,QM3,QM4,AA)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 16  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!     ELEMENTS:                                                    *
!                                                                  *
!       N1       (k1)   N1'         N2       (k2)    N2'        +- *
!     (j  Q J ::W(11)::j  Q'J') * (j  Q J ::W(22):: j  Q'J')    -+ *
!       1  1 1          1  1 1      2  2 2           2  2 2     ++ *
!                                                               -- *
!     SUBROUTINE CALLED: WJ1                                       *
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
      USE CONS_C,          ONLY: ZERO, TENTH, EPS
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K1, K2
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(OUT)              :: AA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4, IQMM12, IQMM34
      REAL(DOUBLE) :: A1, W
!-----------------------------------------------
      AA=ZERO
      IF(IK1(3) > 9) THEN
        IF(IK1(4) > 2) CALL MES(35)
        IF(ID1(4) > 2) CALL MES(35)
      ENDIF
      IF(IK2(3) > 9) THEN
        IF(IK2(4) > 2) CALL MES(35)
        IF(ID2(4) > 2) CALL MES(35)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM12=IQMM1+IQMM2
      IF(IK1(4) /= (ID1(4)+IQMM12))RETURN
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM4=QM4+QM4+TENTH*QM4
      IQMM34=IQMM3+IQMM4
      IF(IK2(4) /= (ID2(4)+IQMM34))RETURN
      CALL WJ1(IK1,BK1,ID1,BD1,K1,QM1,QM2,A1)
      IF(DABS(A1) < EPS)RETURN
      CALL WJ1(IK2,BK2,ID2,BD2,K2,QM3,QM4,W)
      IF(DABS(W) < EPS)RETURN
      AA=A1*W
      RETURN
      END SUBROUTINE GG1122
