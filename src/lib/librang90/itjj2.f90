!*******************************************************************
!                                                                  *
      INTEGER FUNCTION ITJJ2(IK,ID,KG1,BK,BD,IBT,BT,ITP,ITG,IQ)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 07  --------------   *
!                                                                  *
!     FUNCTION CALLED: ITTK                                        *
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
      USE CONS_C,          ONLY: HALF
      USE ribojj_C
      USE ribojj9_C
      USE ribojj11_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)                :: KG1, IQ
      INTEGER,      INTENT(OUT)               :: ITP, ITG
      INTEGER,      INTENT(IN),  DIMENSION(7) :: IK, ID
      INTEGER,      INTENT(OUT), DIMENSION(7) :: IBT
      REAL(DOUBLE), INTENT(IN),  DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT), DIMENSION(3) :: BT
!      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITK, ITD, ITP1, ITG1
!-----------------------------------------------
      ITJJ2=0
      IF(ID(3) > 37) RETURN
      IF(ITTK(ID(6),IK(6),KG1) == 0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3) < 9) THEN
        ITP1=IMPNJJ(ITK)
        ITP=IMPTJJ(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNJJ(ITK)
        ITG=IMGTJJ(ITD)
      ELSEIF(ID(3) == 9) THEN
        IF(ITK > 300) THEN
          IF(ITD < 300) CALL MES(52)
          IF(ID(4) > 2) CALL MES(12)
          IF(IK(4) > 2) CALL MES(12)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPNJJ9(ITK)
          ITP=IMPTJJ9(ITD)
          IF(ITP1 /= ITP)RETURN
          ITG1=IMGNJJ9(ITK)
          ITG=IMGTJJ9(ITD)
        ELSE
          PRINT*, "ERROR in ITJJ2"
          STOP
        END IF
      ELSE
        IF(ID(4) > 2) CALL MES(12)
        IF(IK(4) > 2) CALL MES(12)
        ITP1=IMPNJJ11(ITK)
        ITP=IMPTJJ11(ITD)
        IF(ITP1 /= ITP)RETURN
        ITG1=IMGNJJ11(ITK)
        ITG=IMGTJJ11(ITD)
      ENDIF
      IF(ITG1 /= ITG)RETURN
      ITJJ2=1
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      RETURN
      END FUNCTION ITJJ2
