

!***********************************************************************
!                                                                      *
      SUBROUTINE FRMRWF(INDEX, NSUBS, FILNAM)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, OPENFL.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE DEBUG_C
      USE GRID_C
      USE LEFT_C
      USE ORB_C, ONLY: E, NP, NAK, NH
      USE WAVE_C, ONLY: PZ
      USE WHFROM_C, ONLY: SOURCE
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE intrpq_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NSUBS
      CHARACTER  :: FILNAM*(*)
      INTEGER , INTENT(IN) :: INDEX(NNNW)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, IOS, NPY, NAKY, MY, J, LOC, I, LENTH
      REAL(DOUBLE) :: EY, DNORM
      REAL(DOUBLE), DIMENSION(:), pointer :: PA, QA, RA
      LOGICAL :: FOUND
      CHARACTER :: G92RWF*6
!-----------------------------------------------
!
!

      !FORM = 'UNFORMATTED'
      !STATUS = 'OLD'

      CALL OPENFL (23, FILNAM, 'UNFORMATTED', 'OLD', IERR)
      IF (IERR == 1) THEN
         WRITE (ISTDE, *) 'Error openning file "', FILNAM(1:LEN_TRIM(FILNAM)), &
            '"'
         CLOSE(23)
         STOP
      ENDIF
!
!   Check the file; if not as expected, try again
!
      READ (23, IOSTAT=IOS) G92RWF
      IF (IOS/=0 .OR. G92RWF/='G92RWF') THEN
         WRITE (ISTDE, *) 'This is not a Radial WaveFunction File;'
         CLOSE(23)
         STOP
      ENDIF
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      IF (LDBPR(3)) WRITE (99, 300)
    2 CONTINUE
      FOUND = .FALSE.
      READ (23, IOSTAT=IOS) NPY, NAKY, EY, MY
      IF (IOS == 0) THEN
         DO J = 1, NSUBS
            LOC = INDEX(J)
            IF (.NOT.(.NOT.SET(LOC) .AND. NP(LOC)==NPY .AND. NAK(LOC)==NAKY)) &
               CYCLE
            FOUND = .TRUE.
            E(LOC) = EY
            CALL ALLOC (PA, MY, 'PA', 'FRMFRW')
            CALL ALLOC (QA, MY, 'QA', 'FRMFRW')
            CALL ALLOC (RA, MY, 'RA', 'FRMFRW')
            READ (23) PZ(LOC), (PA(I),I=1,MY), (QA(I),I=1,MY)
            READ (23) (RA(I),I=1,MY)
            CALL INTRPQ (PA, QA, MY, RA, LOC, DNORM)
            IF (LDBPR(3)) WRITE (99, 301) NP(LOC), NH(LOC), E(LOC), DNORM
            CALL DALLOC (PA, 'PA', 'FRMFRW')
            CALL DALLOC (QA, 'QA', 'FRMFRW')
            CALL DALLOC (RA, 'RA', 'FRMFRW')
            LENTH = LEN_TRIM(FILNAM)
            SET(LOC) = .TRUE.
            SOURCE(LOC) = FILNAM(1:3)
            GO TO 2
         END DO
         IF (.NOT.FOUND) THEN
            READ (23)
            READ (23)
            GO TO 2
         ENDIF
      ENDIF
      IF (LDBPR(3)) WRITE (99, *) ' orbitals renormalised;'
!
      CLOSE(23)
!
      RETURN
!
  300 FORMAT(/,'From SUBROUTINE FRMRWF:'/,' Orbital',8X,'Eigenvalue',19X,'Norm'&
         )
  301 FORMAT(2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
      RETURN
!
      END SUBROUTINE FRMRWF
