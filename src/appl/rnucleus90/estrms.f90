!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION ESTRMS (APARM, CPARM) 
!                                                                      *
!   Determines the root mean square radius for a Fermi nucleus given   *
!   the parameters `c' (CPARM) and `a' (APARM). We use the formalism   *
!   developed in F. A. Parpia and A. K. Mohanty ``Relativistic basis   *
!   set calculations for atoms with Fermi nuclei'' Phys Rev A (1992)   *
!   in press.                                                          *
!                                                                      *
!   Call(s) to: SKFUN.                                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:31:28   1/ 3/07  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE skfun_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: APARM 
      REAL(DOUBLE) , INTENT(IN) :: CPARM 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: PI, SQTBF, ABC, PABC, CBAM, DNUMER, DDENOM 
!-----------------------------------------------
!
      PI = 4.0D00*ATAN(1.0D00) 
      SQTBF = SQRT(3.0D00/5.0D00) 
!
      ABC = APARM/CPARM 
      PABC = PI*ABC 
      CBAM = -CPARM/APARM 
      DNUMER = 1.0D00 + (10.0D00/3.0D00)*PABC**2 + (7.0D00/3.0D00)*PABC**4 - &
         120.0D00*ABC**5*SKFUN(5,CBAM) 
      DDENOM = 1.0D00 + PABC**2 - 6.0D00*ABC**3*SKFUN(3,CBAM) 
      ESTRMS = CPARM*SQTBF*SQRT(DNUMER/DDENOM) 
!
      RETURN  
      END FUNCTION ESTRMS 
