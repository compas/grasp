!***********************************************************************
!                                                                      *
      SUBROUTINE ESTIM(J) 
!                                                                      *
!   This  subprogram implements  Part 1 of Algorithm 7.1 of C Froese   *
!   Fischer, Comput Phys Rep, 3 (1986) 320-321.                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE def_C
      USE grid_C
      USE orb_C
      USE scf_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NPJ, NAKABS 
      REAL(DOUBLE) :: ALPHA, CSQ, FNREL, FKABS, FKAP2, ZALPHA, GAMMA, EBYM 
!-----------------------------------------------
!
!   Initializations
!
      ALPHA = 1.0D00/C 
      CSQ = C*C 
      NPJ = NP(J) 
      NAKABS = ABS(NAK(J)) 
      FNREL = DBLE(NPJ - NAKABS) 
      FKABS = DBLE(NAKABS) 
      FKAP2 = FKABS*FKABS 
!
!   Set ZINF, the asymptotic charge seen by the electron
!
!     ZINF = DBLE (IONCTY+1)
!
!   Changed on 07/06/93 by WPW
!
      ZINF = Z + DBLE((-NELEC) + 1) 
!
!   Set the lower bound
!
      ZALPHA = ZINF*ALPHA 
      IF (ZALPHA < FKABS) THEN 
         GAMMA = SQRT(FKAP2 - ZALPHA*ZALPHA) 
         EBYM = 1.0D00/SQRT(1.0D00 + (ZALPHA/(GAMMA + FNREL + 0.5D00))**2) 
         EPSMIN = (1.0D00 - EBYM)*CSQ 
      ELSE 
         EPSMIN = 0.25D00*CSQ/DBLE(NPJ*NPJ) 
      ENDIF 
      EMIN = EPSMIN 
!
!   Set the upper bound
!
      ZALPHA = Z*ALPHA 
      IF (ZALPHA < FKABS) THEN 
         GAMMA = SQRT(FKAP2 - ZALPHA*ZALPHA) 
         EBYM = 1.0D00/SQRT(1.0D00 + (ZALPHA/(GAMMA + FNREL - 0.5D00))**2) 
         EPSMAX = (1.0D00 - EBYM)*CSQ 
      ELSE 
         EPSMAX = CSQ + CSQ 
      ENDIF 
      EMAX = EPSMAX 
!
      RETURN  
      END SUBROUTINE ESTIM 
