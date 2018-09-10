!***********************************************************************
!                                                                      *
      SUBROUTINE CONSIS(J) 
!                                                                      *
!   This routine computes the weighted self-consistency of orbital J   *
!                                                                      *
!   Written by Farid A Parpia, at OXFORD    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE grid_C
      USE int_C
      USE scf_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTP, I 
      REAL(DOUBLE) :: SCMEA, DELTAO 
!-----------------------------------------------
!
!
      SCMEA = 0.0D00 
      MTP = MIN(MTP0,MF(J)) 
      DO I = 1, MTP 
         DELTAO = ABS(P(I)-PF(I,J)) + ABS(Q(I)-QF(I,J)) 
         SCMEA = DMAX1(DELTAO,SCMEA) 
      END DO 
      SCNSTY(J) = SCMEA*SQRT(UCF(J)) 
!
      RETURN  
      END SUBROUTINE CONSIS 
