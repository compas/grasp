

!***********************************************************************
!                                                                      *
      SUBROUTINE FRMTFP(INDEX, NSUBS)
!                                                                      *
!   This  subroutine  is  used  to  produce  estimates  of  the wave   *
!   functions by use of the Thomas-Fermi approximation to the direct   *
!   potential. SUBROUTINE  SOLVH  is used  to obtain the radial wave   *
!   functions.                                                         *
!                                                                      *
!   Call(s) to: [ERWF]: SOLVH.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE LEFT_C, ONLY: SET
      USE ORB_C, ONLY: NP, NH
      USE WAVE_C
      USE WHFROM_C, ONLY: SOURCE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE solvh_I
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NNNP = 590
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NSUBS
      INTEGER , INTENT(IN) :: INDEX(NNNW)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NNN1 = NNNP + 10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, LOC
      LOGICAL :: FAIL
!-----------------------------------------------
!
!
      DO J = 1, NSUBS
!
         LOC = INDEX(J)
!
         IF (SET(LOC)) CYCLE
!
!   Estimate the leading coefficient of the series expansion for
!   the orbital
!
         PZ(LOC) = 10.0D00
!
!   Calculate radial wavefunctions
!
         CALL SOLVH (LOC, FAIL)
!
!   Message if SOLVH did not converge; reset SET(LOC) otherwise
!
         IF (FAIL) THEN
            WRITE (*, 300) NP(LOC), NH(LOC)
         ELSE
            SET(LOC) = .TRUE.
               !SOURCE(LOC) = 'Thomas-Fermi estimate'
            SOURCE(LOC) = 'T-F'
         ENDIF
!
      END DO
!
      RETURN
!
  300 FORMAT(/,'TFWAVE: Unable to compute radial'/,/,' wavefunction for ',I2,A2&
         ,' subshell;')
      RETURN
!
      END SUBROUTINE FRMTFP
