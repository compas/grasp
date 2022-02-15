

!***********************************************************************
!                                                                      *
      SUBROUTINE WRTRWF
!                                                                      *
!   Open, write a header and all subshell radial wavefunctions, and    *
!   close the  .rwf  file.                                             *
!                                                                      *
!   Call(s) to: [LIB92]: DALLOC, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE GRID_C
      USE ORB_C
      USE WAVE_C, ONLY: PZ, PF, QF, MF
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEWUNIT, IERR, J, MFJ, I
      LOGICAL :: IMOPENED
      CHARACTER :: FILNAM*128
!-----------------------------------------------
!
!

      FILNAM = 'rwfn.inp'

      DO NEWUNIT = 23, 99                        ! 23 is a historical value
         INQUIRE(UNIT=NEWUNIT, OPENED=IMOPENED)
         IF (IMOPENED) CYCLE
         EXIT                                    ! should be the normal exit point
      END DO

      IF (NEWUNIT == 100) THEN
         WRITE (ISTDE, *) 'All unit numbers from 23 to 99 are BUSY!'
         STOP
      ENDIF

      CALL OPENFL (NEWUNIT, FILNAM, 'UNFORMATTED', 'NEW', IERR)
      IF (IERR == 1) THEN
         WRITE (ISTDE, *) 'Error when opening "', FILNAM(1:LEN_TRIM(FILNAM)), &
            '"'

         STOP
      ENDIF
!
!   Write the file header
!
      WRITE (NEWUNIT) 'G92RWF'
!
!   Write out the radial wavefunctions
!
      DO J = 1, NW
         MFJ = MF(J)
         WRITE (NEWUNIT) NP(J), NAK(J), E(J), MFJ
         WRITE (NEWUNIT) PZ(J), (PF(I,J),I=1,MFJ), (QF(I,J),I=1,MFJ)
         WRITE (NEWUNIT) (R(I),I=1,MFJ)
      END DO

      CLOSE(NEWUNIT)
!
!   Deallocate the storage for the radial wavefunctions
!
      CALL DALLOC (PF, 'PF', 'WRTRWF')
      CALL DALLOC (QF, 'QF', 'WQRTRWF')
!
      RETURN
      END SUBROUTINE WRTRWF
