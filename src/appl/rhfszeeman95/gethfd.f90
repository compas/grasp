!************************************************************************
!*                                                                      *
       SUBROUTINE GETHFD(NAME)
!*                                                                      *
!*   Interactively determines the data governing the HFS problem.       *
!*                                                                      *
!*   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
!*               [RCI92]: SETISO, SETRWF.                               *
!*                                                                      *
!*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!*                                                                      *
!*   Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03 1/3/07*  
!*   Modified by Charlotte Froese Fischer                               *
!*                     Gediminas Gaigalas  11/01/17                     *
!*   Modified by Wenxian Li   12/28/18                                  *
!************************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW
      USE decide_C
      USE def_C
      USE defaulT_C
      USE foparm_C
      USE grid_C,          ONLY: rnt, h, n, hp
      USE npar_C
      USE orb_C
      USE wave_C
      USE wfac_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s 
!-----------------------------------------------
      USE getyn_I 
      USE setiso_I 
      USE setqic_I 
      USE radgrd_I 
      USE nucpot_I 
      USE setrwfa_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s 
!-----------------------------------------------
      CHARACTER(LEN=24), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s 
!-----------------------------------------------
      LOGICAL            :: YES 
      CHARACTER(LEN=128) :: ISOFILE 
!-----------------------------------------------
!
!   Open, check, load data from, and close the  .iso  file
!
      ISOFILE = 'isodata'
      CALL SETISO (ISOFILE)
!
!   Determine the physical effects specifications
!
      IF (NDEF /= 0) THEN
         WRITE (6, *) 'The physical speed of light in'
         WRITE (6, *) ' atomic units is',CVAC,';'
         WRITE (6, *) ' revise this value?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (6, *) 'Enter the revised value:'
            READ (5, *)C
         ELSE
            C = CVAC
         ENDIF
      ELSE
         C = CVAC
      ENDIF
!
      IF (NDEF /= 0) THEN
!
         WRITE (6, *) 'Treat contributions of some CSFs'
         WRITE (6, *) ' as first-order perturbations?'
         YES = GETYN ()
         IF (YES) THEN
            LFORDR = .TRUE.
            WRITE (6, *) 'The contribution of CSFs'
            WRITE (6, *) ' 1 -- ICCUT will be treated'
            WRITE (6, *) ' variationally; the remainder'
            WRITE (6, *) ' perturbatively; enter ICCUT:'
            READ (5, *) ICCUT
         ELSE
            LFORDR = .FALSE.
            ICCUT = 0
         ENDIF
      ELSE
         LFORDR = .FALSE.
         ICCUT = 0
      ENDIF
!
!   Determine the parameters controlling the radial grid
!
      IF (NPARM == 0) THEN
         RNT = EXP(-65.0D00/16.0D00)/Z
         H = 0.5D00**4
         N = MIN (220,NNNP)
      ELSE
         RNT = 2.0D-06/Z
         H = 5.0D-02
         N = NNNP
      ENDIF
      HP = 0.0D00
      IF (NDEF /= 0) THEN
         WRITE (6, *) 'The default radial grid parameters'
         WRITE (6, *) ' for this case are:'
         WRITE (6, *) ' RNT = ',RNT,';'
         WRITE (6, *) ' H = ',H,';'
         WRITE (6, *) ' HP = ',HP,';'
         WRITE (6, *) ' N = ',N,';'
         WRITE (6, *) ' revise these values?'
         YES = GETYN ()
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
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD
!
!   Generate $- r \times V_nuc (r)$
!
      CALL NUCPOT
!
!   Load the radial wavefunctions
!
!c      CALL SETRWFA(NAME)

      CALL SETRWFA(TRIM(NAME)//'.w')

!
      RETURN
      END SUBROUTINE GETHFD
