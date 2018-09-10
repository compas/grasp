!*******************************************************************
!                                                                  *
      SUBROUTINE PERKO1(JA,BK,IK,BD,ID)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 22  -------------  *
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (FOR ONE SHELL)    *
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
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: HALF
      USE m_C,             ONLY: JLIST, NQ1, NQ2, JJQ1, JJQ2
      USE orb_C,           ONLY: NP, NAK
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nmtejj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)                :: JA
      INTEGER,      INTENT(OUT), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(OUT), DIMENSION(3) :: BK, BD
!      DIMENSION BK(3),IK(7),BD(3),ID(7)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJ
!-----------------------------------------------
      IJ=JLIST(JA)
      IK(2)=NP(IJ)
      ID(2)=IK(2)
      IK(3)=(IABS(NAK(IJ))*2)-1
      ID(3)=IK(3)
      IK(4)=NQ1(IJ)
      ID(4)=NQ2(IJ)
      IK(5)=(IK(3)+NAK(IJ)/IABS(NAK(IJ)))/2
      ID(5)=IK(5)
      IK(6)=JJQ1(3,IJ)-1
      ID(6)=JJQ2(3,IJ)-1
      IK(7)=IABS(NAK(IJ))-JJQ1(1,IJ)
      ID(7)=IABS(NAK(IJ))-JJQ2(1,IJ)
      BK(1)=HALF*DBLE(IK(7))
      BD(1)=HALF*DBLE(ID(7))
      BK(2)=HALF*DBLE(IK(6))
      BD(2)=HALF*DBLE(ID(6))
      BK(3)=-HALF*DBLE(IABS(NAK(IJ))-IK(4))
      BD(3)=-HALF*DBLE(IABS(NAK(IJ))-ID(4))
      IK(1)=NMTEJJ(IK(7),IK(6),IK(3),ID(4),IK(4))
      ID(1)=NMTEJJ(ID(7),ID(6),ID(3),ID(4),IK(4))
      RETURN
      END SUBROUTINE PERKO1
