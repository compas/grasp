 
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION OUTBND (ETRY) 
!                                                                      *
!   This  subprogram determines whether the trial eigenvalue etry is   *
!   within the bounds (EPSMIN,EPSMAX)                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE SCF_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: ETRY 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
      IF (ETRY>EPSMIN .AND. ETRY<EPSMAX) THEN 
         OUTBND = .FALSE. 
      ELSE 
         OUTBND = .TRUE. 
      ENDIF 
!
      RETURN  
      END FUNCTION OUTBND 
