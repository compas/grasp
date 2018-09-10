!*******************************************************************
!                                                                  *
      SUBROUTINE GG1112(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,K1,  &
                        QM1,QM2,QM3,QM4,WW)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 15  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!     ELEMENT:                                                     *
!                                                                  *
!       N1       (k1) (j1)(j2)  N1'       N2     (j2)   N2'     +- *
!     (j Q J ::[W(11)*A(1)]  ::j  Q'J')*(j Q J ::A(2)::j Q'J')  -+ *
!       1 1 1                   1  1 1    2 2 2         2 2 2   ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ,WAP1                         *
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
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C,          ONLY: ZERO, TENTH, HALF, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE c0t5s_I
      USE rmeajj_I
      USE wap1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK1, IK2, ID1, ID2
      REAL(DOUBLE), INTENT(IN)               :: QM1, QM2, QM3, QM4
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK1, BK2, BD1, BD2
      REAL(DOUBLE), INTENT(OUT)              :: WW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4, IQMM23, KK1
      REAL(DOUBLE) :: A1, S, BK, WA
!-----------------------------------------------
      WW=ZERO
      IF(IK1(3) > 9) THEN
        IF(IK1(4) > 2) CALL MES(34)
        IF(ID1(4) > 2) CALL MES(34)
      ENDIF
      IF(IK2(3) > 9) THEN
        IF(IK2(4) > 2) CALL MES(34)
        IF(ID2(4) > 2) CALL MES(34)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM23=IQMM1+IQMM2+IQMM3
      IF(IK1(4) /= (ID1(4)+IQMM23))RETURN
      IQMM4=QM4+QM4+TENTH*QM4
      IF(IK2(4) /= (ID2(4)+IQMM4))RETURN
      KK1=K1*2
      CALL C0T5S(BD2(1),BD2(3),QM4,BK2(1),BK2(3),A1)
      IF(ABS(A1) < EPS)RETURN
!GG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
      IF(ABS(S) < EPS)RETURN
      BK=HALF*DBLE(IK2(3))
      CALL WAP1(IK1,BK1,ID1,BD1,K1,BK,QM1,QM2,QM3,WA)
      IF(ABS(WA) < EPS)RETURN
      WW=-A1*WA*S/SQRT(DBLE(IK2(7)+1))
      RETURN
      END SUBROUTINE GG1112
