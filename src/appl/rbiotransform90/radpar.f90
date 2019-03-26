!***********************************************************************
!                                                                      *
      SUBROUTINE RADPAR
!                                                                      *
!   This subroutine sets the parameters controlling the radial grid    *
!                                                                      *
!   Last revision: June 1996                                           *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE def_C,           ONLY: c, cvac, z, accy
      USE default_C
      USE grid_C
      USE npar_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL   :: YES
      CHARACTER :: ANSWER
!-----------------------------------------------
!
      IF (NPARM == 0) THEN
         RNT = EXP((-65.0D00/16.0D00))/Z
         H = 0.5D00**4
         N = MIN(220,NNNP)
      ELSE
!cff     .. Should be Z-dependent
         RNT = 2.0D-06/Z
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D00
      IF (NDEF /= 0) THEN
         WRITE (6, *) 'The default radial grid parameters'
         WRITE (6, *) ' for this case are:'
         WRITE (6, *) ' RNT = ', RNT, ';'
         WRITE (6, *) ' H = ', H, ';'
         WRITE (6, *) ' HP = ', HP, ';'
         WRITE (6, *) ' N = ', N, ';'
         WRITE (6, *) ' revise these values?'
         YES = GETYN()
         IF (YES) THEN
            WRITE (6, *) 'Enter RNT:'
            READ (5, *) RNT
            WRITE (6, *) 'Enter H:'
            READ (5, *) H
            WRITE (6, *) 'Enter HP:'
            READ (5, *) HP
            WRITE (6, *) 'Enter N:'
            READ (5, *) N
         ENDIF
      ENDIF
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!

      IF (NDEF /= 0) THEN
         WRITE (6, *) 'The physical speed of light in'
         WRITE (6, *) ' atomic units is', CVAC, ';'
         WRITE (6, *) ' revise this value?'
         YES = GETYN()
         IF (YES) THEN
            WRITE (6, *) 'Enter the revised value:'
            READ (5, *) C
         ELSE
            C = CVAC
         ENDIF
      ELSE
         C = CVAC
      ENDIF

      RETURN
      END SUBROUTINE RADPAR
