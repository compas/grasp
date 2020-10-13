

!***********************************************************************
!                                                                      *
SUBROUTINE FRMHYD(INDEX, NSUBS, MODIFY)
!                                                                      *
!   This  subroutine  is  used  to  produce  estimates  of  the wave   *
!   functions from the hydrogenic approximation. The screening cons-   *
!   tant is taken to be 0.4 (Changed, see comments below,XHH)          *
!                                                                      *
!   Call(s) to: [ERWF]: DCBSRW.                                        *
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
      USE DEF_C
      USE LEFT_C
      USE ORB_C, ONLY: E, NP, NAK, NH
      USE WAVE_C, ONLY: PF, QF, PZ, MF
      USE WHFROM_C, ONLY: SOURCE
      USE HYDPAR_C, ONLY: SIGMA
      USE IOUNIT_C
      USE DEFAULT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dcbsrw_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NSUBS
      LOGICAL , INTENT(IN) :: MODIFY
      INTEGER , INTENT(IN) :: INDEX(NNNW)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      !INTEGER, PARAMETER :: NNNP = 590
      INTEGER, PARAMETER :: NNN1 = NNNP + 10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, LOC
      REAL(DOUBLE) :: ZEFF
      LOGICAL :: YES, GETYN
!-----------------------------------------------
!
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) '***** Screening parameters ******'
      DO J = 1, NSUBS
         LOC = INDEX(J)
         IF (SET(LOC)) CYCLE

         WRITE (ISTDE, '(I2,A2,F10.2)') NP(LOC), NH(LOC), SIGMA(LOC)
         IF (MODIFY) THEN
            WRITE (ISTDE, *) 'Input new value >'
            READ (5, *) SIGMA(LOC)
         ENDIF

         ZEFF = Z - SIGMA(LOC)

!           ...Calculate radial wavefunctions
         CALL DCBSRW (NP(LOC), NAK(LOC), ZEFF, E(LOC), PZ(LOC), PF(:,LOC), &
            QF(:,LOC), MF(LOC))
         SET(LOC) = .TRUE.
            !SOURCE(LOC) = 'Screened hydrogenic estimate'
         SOURCE(LOC) = 'Hyd'
      END DO
      RETURN
      END SUBROUTINE FRMHYD
