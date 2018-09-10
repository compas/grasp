!***********************************************************************
!                                                                      *
      SUBROUTINE GETCPR(RRMS, APARM, CPARM) 
!                                                                      *
!   Determines the parameter `c' (CPARM) for a Fermi nucleus,  given   *
!   the root mean square radius (RRMS) and the parameter `a' (APARM).  *
!   We use the formalism developed in F. A. Parpia and A. K. Mohanty   *
!   ``Relativistic basis set  calculations for atoms with  Fermi nu-   *
!   clei'' Phys Rev A (1992) in press.                                 *
!                                                                      *
!   Call(s) to: ESTRMS.                                                *
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
      USE estrms_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: RRMS 
      REAL(DOUBLE)  :: APARM 
      REAL(DOUBLE) , INTENT(OUT) :: CPARM 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: ACCY, CPMIN, CPMAX, CPTRY, RMSTRY 
!-----------------------------------------------
!
!xhb
!xh Accuracy parameter
!xhe
      ACCY = 1.0D-12 
!
!   Bracket CPARM with a lower and upper limit
!
!   Lower limit
!
      CPMIN = RRMS 
      CPMIN = 0.5D00*CPMIN 
      DO WHILE(ESTRMS(APARM,CPMIN) > RRMS) 
         CPMIN = 0.5D00*CPMIN 
      END DO 
!
!   Upper limit
!
      CPMAX = RRMS 
      CPMAX = 2.0D00*CPMAX 
      DO WHILE(ESTRMS(APARM,CPMAX) < RRMS) 
         CPMAX = 2.0D00*CPMAX 
      END DO 
!
!   Find CPARM by the method of bisection
!
      CPTRY = 0.5D00*(CPMAX + CPMIN) 
!
      RMSTRY = ESTRMS(APARM,CPTRY) 
!
      IF (RMSTRY > RRMS) THEN 
         CPMAX = CPTRY 
      ELSE 
         CPMIN = CPTRY 
      ENDIF 
      DO WHILE((CPMAX - CPMIN)/(CPMAX + CPMIN)>ACCY .AND. ABS(RMSTRY-RRMS)/RRMS&
         >ACCY) 
         CPTRY = 0.5D00*(CPMAX + CPMIN) 
!
         RMSTRY = ESTRMS(APARM,CPTRY) 
!
         IF (RMSTRY > RRMS) THEN 
            CPMAX = CPTRY 
         ELSE 
            CPMIN = CPTRY 
         ENDIF 
!
      END DO 
!
      CPARM = CPTRY 
!
      RETURN  
      END SUBROUTINE GETCPR 
