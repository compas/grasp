!************************************************************************
!*                                                                      *
       SUBROUTINE STRSUM
!*                                                                      *
!*   Generates the first part of <name>.ch or <name>.h on stream 29     *
!*                                                                      *
!*   Call(s) to: ENGOUTH                                                *
!*                                                                      *
!*   Written by P. Jonsson                 Last revision: 20 Oct 1999   *
!*                                                                      *
!*   Translated by Wenxian Li F77 to F90 12/28/18                       *
!************************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE eigv_C
      USE nsmdat_C, ONLY: SQN, DMOMNM, QMOMB, HFSI=>SQN,   &
                          HFSD=>DMOMNM, HFSQ=>QMOMB
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
      WRITE (29,302) HFSI
      WRITE (29,303) HFSD
      WRITE (29,304) HFSQ
!
!    Write the list of eigenpair indices
!
      CALL ENGOUTH (EAV, EVAL, IATJPO, IASPAR, IVEC, NVEC, 3)
!
!
  302 FORMAT ('Nuclear spin                        ',1PD22.15,' au')
  303 FORMAT ('Nuclear magnetic dipole moment      ',1PD22.15,' n.m.')
  304 FORMAT ('Nuclear electric quadrupole moment  ',1PD22.15,' barns')
      RETURN
!
      END SUBROUTINE STRSUM
