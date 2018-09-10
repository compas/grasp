!*******************************************************************
!                                                                  *
      SUBROUTINE ONEPARTICLEJJ2(NS,KA,JA,JB,COEFF)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1 +- 1   *
!                                                  N'2 = N2 -+ 1   *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
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
      USE m_C,             ONLY: NQ1, JLIST
      USE orb_C,           ONLY: NAK
      USE trk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recop00_I
      USE recop2_I
      USE c0t5s_I
      USE rmeajj_I
      USE wj1_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NS,KA,JA,JB
      REAL(DOUBLE), INTENT(OUT) :: COEFF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: I,IA,IB,IAT,IJ,IFAZ,IQMM1,IQMM2,JIBKS1,JIB, &
                      KS1,KS2
      REAL(DOUBLE) :: REC,A1,A2,A3,S1,S2,QM1,QM2
!-----------------------------------------------
!
!     THE CASE 12   + -
!
      COEFF=ZERO
      IA=MIN0(JA,JB)
      IB=MAX0(JA,JB)
      CALL RECOP00(NS,IA,IB,KA,IAT)
      IF(IAT == 0)RETURN
      IJ=JLIST(IA)
      KS1=(IABS(NAK(IJ))*2)-1
      IJ=JLIST(IB)
      KS2=(IABS(NAK(IJ))*2)-1
      CALL RECOP2(NS,IA,IB,KS1,KS2,KA,0,IAT,REC)
      IF(IAT == 0)RETURN
      CALL PERKO2(JA,JB,JA,JA,2)
      QM1=HALF
      QM2=-HALF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4) /= (ID1(4)+IQMM1)) RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4) /= (ID2(4)+IQMM2)) RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A2)
      IF(DABS(A2) < EPS) RETURN
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),A3)
      IF(DABS(A3) < EPS) RETURN
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S1)
      IF(DABS(S1) < EPS) RETURN
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S2)
      IF(DABS(S2) < EPS) RETURN
      A1=S1*S2*A2*A3
      CALL RECOP2(NS,IA,IB,KS1,KS2,KA,1,IAT,REC)
      COEFF=A1*REC/DSQRT(DBLE((2*KA+1)*(IK1(7)+1)*(IK2(7)+1)))
      JIB=IB-1
      IFAZ=0
      DO  I=IA,JIB
        IJ=JLIST(I)
        IFAZ=IFAZ+NQ1(IJ)
      END DO
      IFAZ=IFAZ+1
      IF(MOD(IFAZ,2) /= 0)COEFF=-COEFF
      IF(IA.NE.JA) THEN
        IFAZ = IK1(3)+IK2(3)-2*KA+2
        IF(MOD(IFAZ,4) /= 0)COEFF=-COEFF
      ENDIF
      COEFF=-COEFF*SQRT(DBLE(ID1(3)+1))
      RETURN
      END SUBROUTINE ONEPARTICLEJJ2
