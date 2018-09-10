!*******************************************************************
!                                                                  *
      SUBROUTINE GG1222(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,K1, &
                        QM1,QM2,QM3,QM4,WW)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 17  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                                          N1      (j1)   N1'      *
!     ELEMENTS:                          (j  Q J ::A(1)::j Q'J')*  *
!                                          1  1 1         1 1 1    *
!        N2        (j2)    (k1) (j1)  N2'                       +- *
!     *(j  Q J ::[ A(2) * W(22) ]  ::j  Q'J')                   -+ *
!        2  2 2                       2  2 2                    ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ,AWP1                         *
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
      USE CONS_C,          ONLY: ZERO, TENTH, HALF, EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE c0t5s_I
      USE rmeajj_I
      USE awp1_I
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
      INTEGER      :: IQMM1,IQMM2,IQMM3,IQMM4,IQMM34,KK1
      REAL(DOUBLE) :: A1, S, BK, AW
!-----------------------------------------------
      WW=ZERO
      IF(IK1(3) > 9) THEN
        IF(IK1(4) > 2) CALL MES(33)
        IF(ID1(4) > 2) CALL MES(33)
      ENDIF
      IF(IK2(3) > 9) THEN
        IF(IK2(4) > 2) CALL MES(33)
        IF(ID2(4) > 2) CALL MES(33)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM4=QM4+QM4+TENTH*QM4
      IQMM34=IQMM2+IQMM3+IQMM4
      IF(IK2(4) /= (ID2(4)+IQMM34))RETURN
      KK1=K1*2
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(DABS(A1) < EPS)RETURN
!GG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      IF(DABS(S) < EPS)RETURN
      BK=HALF*DBLE(IK1(3))
      CALL AWP1(IK2,BK2,ID2,BD2,K1,BK,QM2,QM3,QM4,AW)
      IF(DABS(AW) < EPS)RETURN
      WW=-A1*AW*S/DSQRT(DBLE(IK1(7)+1))
      RETURN
      END SUBROUTINE GG1222
