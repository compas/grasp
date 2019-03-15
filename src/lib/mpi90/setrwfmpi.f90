!***********************************************************************
      SUBROUTINE SETRWFmpi (NAME)

!   Open, check, load data from and close the  .rwf  file.             *
!
!   Used by rcimpivu, rscfmpi, rcimpi
!                                                                      *
!   Call(s) to: [LIB92]: LODRWFmpi, OPENFL.                            *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   MPI version by Xinghong He            Last revision: 06 Aug 1998   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mpi_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodrwfmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = *) NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOS, IERROR
      CHARACTER (LEN = 6) :: G92RWF
!-----------------------------------------------
      IF (myid .EQ. 0) THEN
         CALL openfl (23, name, 'UNFORMATTED', 'OLD', ierror)
         IF (ierror .EQ. 1) THEN
            WRITE (istde,*) 'Error opening', name(1:LEN_TRIM (name))
            STOP
         ENDIF

!   Check the file; if not as expected, stop.

         READ (23,IOSTAT = IOS) G92RWF
         IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
            WRITE (istde,*) 'This is not a Radial WaveFunction File;'
            CLOSE (23)
            STOP
         ENDIF
      ENDIF

!   Attempt to load the radial wavefunctions

      CALL LODRWFmpi (ierror)

      IF (ierror .NE. 0) THEN
         IF (myid .EQ. 0) THEN
            WRITE (istde,*) 'Radial wavefunctions defined in CSL file' &
                          , ' not found in Radial WaveFunction File'
            CLOSE (23)
         ENDIF
         STOP
      ENDIF

      IF (myid .EQ. 0) CLOSE (23)

      RETURN
      END
