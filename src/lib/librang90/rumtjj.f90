!*******************************************************************
!                                                                  *
      SUBROUTINE RUMTJJ(KNT,JJ,LQ,LV,L)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 16  --------------   *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
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
      USE mtjj_C
      USE mtjj2_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE jthn_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JJ, KNT
      INTEGER, INTENT(OUT) :: LQ, LV, L
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KT, KNTMIN
!-----------------------------------------------

      IF(JJ < 9) THEN
        KT=MT(KNT)
      ELSE IF(JJ == 9) THEN
        IF(KNT > 300) THEN
          KNTMIN=KNT-300
          KT=MT9(KNTMIN)
        ELSE
          PRINT*, "ERROR in RUMTJJ"
          STOP
!GG          CALL RUMT67(KNT,NR,LQ,LS,L)
!GG          RETURN
        ENDIF
      ELSE
        KT=MT11(KNT)
      ENDIF
      LQ=JTHN(KT,3,100)
      LV=JTHN(KT,2,100)
      L=JTHN(KT,1,100)
      RETURN
      END SUBROUTINE RUMTJJ
