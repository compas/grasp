!***********************************************************************
!                                                                      *
      SUBROUTINE STRSUM 
!                                                                      *
!   Generates the first part of  hfs92.sum  (on stream 24 and 29).     *
!                                                                      *
!   Call(s) to: ENGOUTH                                                *
!                                                                      *
!   Written by P. Jonsson                 Last revision: 20 Oct 1999   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE eigv_C
      USE nsmdat_C,        ONLY: SQN, DMOMNM, QMOMB,                  &
                                 HFSI=>SQN, HFSD=>DMOMNM, HFSQ=>QMOMB
      USE prnt_C
      USE syma_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE engouth_I 
      IMPLICIT NONE
!-----------------------------------------------
!
!    Write nuclear data
!
      WRITE (24, 302) HFSI 
      WRITE (24, 303) HFSD 
      WRITE (24, 304) HFSQ 
      WRITE (29, 302) HFSI 
      WRITE (29, 303) HFSD 
      WRITE (29, 304) HFSQ 
!
!    Write the list of eigenpair indices
!
!      CALL ENGOUTH (EAV, EVAL, IATJPO, IASPAR, IVEC, NVEC, 3) 
!
      RETURN  
!
  302 FORMAT('Nuclear spin                        ',1P,D22.15,' au') 
  303 FORMAT('Nuclear magnetic dipole moment      ',1P,D22.15,' n.m.') 
  304 FORMAT('Nuclear electric quadrupole moment  ',1P,D22.15,' barns') 
      RETURN  
!
      END SUBROUTINE STRSUM 
