!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ICHKQ1 (JA, JB)
!                                                                      *
!   This routine is to check the occupation condition for one electron *
!   operator.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!                                                                      *
!   Yu Zou                                Last revision: 8/16/00       *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:20   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEBUG_C
      USE ORB_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: JA
      INTEGER :: JB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I, I_IQA, I_IQB

!-----------------------------------------------
!
!
      ICHKQ1 = 0
      K = 0
      DO I = 1, NW
         I_IQA = IQ(I,JA)
         I_IQB = IQ(I,JB)
         IF (I_IQA == I_IQB) CYCLE
         K = K + 1
         IF (K > 2) RETURN
         IF (IABS(I_IQA - I_IQB) <= 1) CYCLE
         RETURN
      END DO
      IF (K==2 .OR. K==0) ICHKQ1 = 1
      RETURN
      END FUNCTION ICHKQ1
