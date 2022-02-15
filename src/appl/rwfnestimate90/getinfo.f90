

!***********************************************************************
!                                                                      *
      SUBROUTINE GETINF
!                                                                      *
!   Interactively determines data useful for generating estimates of   *
!   the subshell radial wavefunctions.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE IOUNIT_C
      USE DEF_C
      USE DEFAULT_C
      USE GRID_C, ONLY: RNT, H, HP, N
      USE NPAR_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setiso_I
      USE setqic_I
      USE radgrd_I
      USE nucpot_I
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: YES
!-----------------------------------------------
!
!
!
!   Open, check, load data from, and close the  .iso  file
!
      CALL SETISO ('isodata')
!
!   Set defaults
!
      C = CVAC
!
      IF (NPARM == 0) THEN
         RNT = EXP((-65.0D00/16.0D00))/Z
         H = 0.5D00**4
         N = MIN(220,NNNP)
      ELSE
!CFF     .. should be Z-dependent
         RNT = 2.0D-06/Z
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D00
!
      IF (NDEF /= 0) THEN
         WRITE (ISTDE, *) 'Change the default speed of light ', &
            'or radial grid parameters?'
!
         YES = GETYN()
         IF (YES) THEN
!
!   Modify the speed of light
!
            WRITE (ISTDE, *) 'The physical speed of light in ', &
               'atomic units is', CVAC, ';'
            WRITE (ISTDE, *) 'revise this value?'
            YES = GETYN()
            IF (YES) THEN
               WRITE (ISTDE, *) 'Enter the revised value:'
               READ (5, *) C
            ENDIF
!
!   Modify the parameters controlling the radial grid
!
            WRITE (ISTDE, *) 'The default radial grid parameters ', &
               'for this case are:'
            WRITE (ISTDE, *) ' RNT = ', RNT, ';'
            WRITE (ISTDE, *) ' H = ', H, ';'
            WRITE (ISTDE, *) ' HP = ', HP, ';'
            WRITE (ISTDE, *) ' N = ', N, ';'
            WRITE (ISTDE, *) 'revise these values?'
            YES = GETYN()
            IF (YES) THEN
               WRITE (ISTDE, *) 'Enter RNT:'
               READ (5, *) RNT
               WRITE (ISTDE, *) 'Enter H:'
               READ (5, *) H
               WRITE (ISTDE, *) 'Enter HP:'
               READ (5, *) HP
               WRITE (ISTDE, *) 'Enter N:'
               READ (5, *) N
            ENDIF
!
         ENDIF
      ENDIF
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD
!
!   Generate $- r \times V_ (r)$
!
      CALL NUCPOT
!
      RETURN
      END SUBROUTINE GETINF
