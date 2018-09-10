!*******************************************************************
!                                                                  *
      SUBROUTINE GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
      ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,WW)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 19  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                      N1     (j1)  N1'      N2    (j2)   N2'      *
!     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
!                      1 1 1        1 1 1    2 2 2        2 2 2    *
!                                                                  *
!        N3     (j3)  N3'      N4     (j4)  N4'                 +- *
!     *(j Q J::A(3)::j Q'J')*(j Q J::A(4)::j Q'J')              -+ *
!        3 3 3        3 3 3    4 4 4        4 4 4               ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ                              *
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE c0t5s_I
      USE rmeajj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK1,IK2,IK3,IK4,ID1,ID2,ID3,ID4
      REAL(DOUBLE), INTENT(IN)               :: QM1,QM2,QM3,QM4
      REAL(DOUBLE), INTENT(IN), DIMENSION(3) :: BK1,BK2,BK3,BK4,BD1,BD2,BD3,BD4
      REAL(DOUBLE), INTENT(OUT)              :: WW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4
      REAL(DOUBLE) :: A1, C, C1, S, V, Z
!-----------------------------------------------
      WW=ZERO
      IF(IK1(3) > 9) THEN
        IF(IK1(4) > 2) CALL MES(30)
        IF(ID1(4) > 2) CALL MES(30)
      ENDIF
      IF(IK2(3) > 9) THEN
        IF(IK2(4) > 2) CALL MES(30)
        IF(ID2(4) > 2) CALL MES(30)
      ENDIF
      IF(IK3(3) > 9) THEN
        IF(IK3(4) > 2) CALL MES(30)
        IF(ID3(4) > 2) CALL MES(30)
      ENDIF
      IF(IK4(3) > 9) THEN
        IF(IK4(4) > 2) CALL MES(30)
        IF(ID4(4) > 2) CALL MES(30)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4) /= (ID2(4)+IQMM2))RETURN
      IQMM3=QM3+QM3+TENTH*QM3
      IF(IK3(4) /= (ID3(4)+IQMM3))RETURN
      IQMM4=QM4+QM4+TENTH*QM4
      IF(IK4(4) /= (ID4(4)+IQMM4))RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(DABS(A1) < EPS)RETURN
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
      IF(DABS(C1) < EPS)RETURN
      A1=A1*C1
      CALL C0T5S(BD3(1),BD3(3),QM3,BK3(1),BK3(3),C1)
      IF(DABS(C1) < EPS)RETURN
      A1=A1*C1
      CALL C0T5S(BD4(1),BD4(3),QM4,BK4(1),BK4(3),C1)
      IF(DABS(C1) < EPS)RETURN
      A1=A1*C1
!GG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      IF(DABS(S) < EPS)RETURN
!GG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      IF(DABS(C) < EPS)RETURN
!GG      CALL SJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
      CALL RMEAJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
      IF(DABS(V) < EPS)RETURN
!GG      CALL SJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
      CALL RMEAJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
      IF(DABS(Z) < EPS)RETURN
      WW=A1*S*C*V*Z/DSQRT(DBLE((IK1(7)+1)*(IK2(7)+1)*(IK3(7)+1)*(IK4(7)+1)))
      RETURN
      END SUBROUTINE GG1234
