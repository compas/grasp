!***********************************************************************
!                                                                      *
      SUBROUTINE SPEAK(JA, JB, IA1, IB1, IA2, IB2, K, X) 
!                                                                      *
!   Output MCP  coefficients and integral parameters to COMMON block   *
!   /BUFFER/. Also print these if  IBUG1 = 1 .                         *
!                                                                      *
!                                           Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:59:36   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE BUFFER_C,        ONLY: LABEL, COEFF, NBDIM, NVCOEF
      USE DEBUG_C ,        ONLY: IBUG1
      USE ORB_C,           ONLY: NP, NH 
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA 
      INTEGER, INTENT(IN) :: JB 
      INTEGER, INTENT(IN) :: IA1 
      INTEGER, INTENT(IN) :: IB1 
      INTEGER, INTENT(IN) :: IA2 
      INTEGER, INTENT(IN) :: IB2 
      INTEGER, INTENT(IN) :: K 
      REAL(DOUBLE), INTENT(IN) :: X 
!-----------------------------------------------
!
!
      IF (IBUG1 /= 0) WRITE (99, 300) JA, JB, NP(IA1), NH(IA1), NP(IB1), NH(IB1&
         ), NP(IA2), NH(IA2), NP(IB2), NH(IB2), K, X 
!
!   Increment counter
!
      NVCOEF = NVCOEF + 1 
!
!   Ensure that arrays are of adequate size; reallocate if necessary
!
      IF (NVCOEF > NBDIM) CALL ALCBUF (2) 
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
      LABEL(1,NVCOEF) = IA1 
      LABEL(2,NVCOEF) = IB1 
      LABEL(3,NVCOEF) = IA2 
      LABEL(4,NVCOEF) = IB2 
      LABEL(5,NVCOEF) = K 
      COEFF(NVCOEF) = X 
!
      RETURN  
!
  300 FORMAT(2(1X,1I2),4(1X,I2,A2),1X,I2,1X,1P,D19.12) 
      RETURN  
!
      END SUBROUTINE SPEAK 
