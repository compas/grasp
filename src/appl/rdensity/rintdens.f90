!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINTDENS (I, J) 
!                                                                      *
!   The value of RINTDENS is an approximation to:                      *
!                                                                      *
!                                                                      *
!                r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
!                     I     J       I     J                            *
!                r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
!                     I     J       I     J                            *
!   Call(s) to: [SMS92]:  POLINT                                       *
!                                                                      *
!   Written by Per Jonsson             Last revision: 24 Dec 1992      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE grid_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE polint_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L 
      REAL(DOUBLE), DIMENSION(3) :: XA, YA 
      REAL(DOUBLE) :: DENS, DDENS 
!-----------------------------------------------
!
      DO L = 4, 2, -1 
         XA(L-1) = R(L) 
         YA(L-1) = (PF(L,I)*PF(L,J) + QF(L,I)*QF(L,J))/(R(L)*R(L)) 
      END DO 
      CALL POLINT (XA, YA, 3, 0.0D00, DENS, DDENS) 
      RINTDENS = DENS 
!
      RETURN  
      END FUNCTION RINTDENS 
