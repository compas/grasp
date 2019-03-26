!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ICHOP (ISUBSH, ICSF)
!                                                                      *
!   ICHOP is -1 if subshell ISUBSH is empty in CSF  ICSF,  +1 if the   *
!   subshell is full, and 0 if it is open.                             *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:25   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE ORB_C,         ONLY: NKL, NKJ
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ISUBSH
      INTEGER :: ICSF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOCC, IFULL
!-----------------------------------------------
!
!
! cff  Since ICHOP is always called from within a do-loop over the
!      appropriate range, testing seems redundant
!     IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NW)) THEN
!        IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
      IOCC = IQ(ISUBSH,ICSF)
      IFULL = NKJ(ISUBSH) + 1
      IF (IOCC == 0) THEN
         ICHOP = -1
      ELSE IF (IOCC == IFULL) THEN
         ICHOP = 1
      ELSE
         ICHOP = 0
      ENDIF
!        ELSE
!           PRINT *, 'ICHOP: Argument ICSF is out of range.'
!           STOP
!        ENDIF
!     ELSE
!        PRINT *, 'ICHOP: Argument ISUBSH is out of range.'
!        STOP
!     ENDIF
!
      RETURN
      END FUNCTION ICHOP
