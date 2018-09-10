!***********************************************************************
!                                                                      *
      SUBROUTINE SETXUV(J) 
!                                                                      *
!   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
!   subprograms IN and  OUT.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE def_C
      USE grid_C
      USE int_C
      USE pote_C
      USE scf_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, I 
      REAL(DOUBLE) :: DMHH 
!-----------------------------------------------
!
!   Define constants
!
      DMHH = -H*0.5D00 
!
!   Set up arrays XU and XV; since XU(1), XV(1) are never used,
!   set them to some arbitrary value
!
      NM1 = N - 1 
      XU(1) = 0.0D00 
      XV(1) = 0.0D00 
      XU(2:NM1) = DMHH*(XP(3:NM1+1)*RPOR(3:NM1+1)+XP(2:NM1)*RPOR(2:NM1)) + DP(2&
         :NM1) 
      XV(2:NM1) = DMHH*(XQ(3:NM1+1)*RPOR(3:NM1+1)+XQ(2:NM1)*RPOR(2:NM1)) + DQ(2&
         :NM1) 
!
      XU(N) = DMHH*XP(N)*RPOR(N) + DP(N) 
      XV(N) = DMHH*XP(N)*RPOR(N) + DQ(N) 
!
      RETURN  
      END SUBROUTINE SETXUV 
