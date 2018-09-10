!***********************************************************************
!                                                                      *
      SUBROUTINE SETMC 
!                                                                      *
!   This subprogram sets machine-dependent parameters.                 *
!                                                                      *
!   Call(s) to: [LAPACK]: DLAMCH.                                      *
!                                                                      *
!   Written by Farid A Parpia              Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:37   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE DEBUG_C 
      USE DEF_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: DNUM 
      REAL(DOUBLE) :: DLAMCH
!-----------------------------------------------
!
!
!   Set the machine-dependent parameters:
!
!   TENMAX - Maximum size of exponent of 10
!
      TENMAX = DLAMCH('L') 
!
!   EXPMAX - Maximum size of exponent of e
!
      DNUM = DLAMCH('O') 
      EXPMAX = LOG(DNUM) 
!
!   EXPMIN - Minimum size of exponent of e
!
      DNUM = DLAMCH('U') 
      EXPMIN = LOG(DNUM) 
!
!   PRECIS - Machine precision
!
      PRECIS = DLAMCH('E') 
!
!   Debug printout
!
      IF (LDBPG(1)) WRITE (99, 300) TENMAX, EXPMAX, EXPMIN, PRECIS 
!
      RETURN  
!
  300 FORMAT(/,'From SUBROUTINE SETMC:'/,' TENMAX (maximum exponent of 10): ',&
         F5.0,/,' EXPMAX (maximum exponent of e): ',1P,1D19.12,/,&
         ' EXPMIN (minimum exponent of e): ',1D19.12,/,&
         ' PRECIS (machine precision): ',1D19.12) 
      RETURN  
!
      END SUBROUTINE SETMC 
