!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ITJPOR (ICSF)
!                                                                      *
!   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE STAT_C,        ONLY: JCUPAR
      USE IOUNIT_C,      ONLY: ISTDE
      USE orb_C,         ONLY: NCFR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ICSF
!-----------------------------------------------
!
      IF (ICSF>=1 .AND. ICSF<=NCFR) THEN
        itjpor = jcupar(NNNW,icsf)
         IF (ITJPOR > 127) ITJPOR = 256 - ITJPOR
         ITJPOR = IABS (ITJPOR)
      ELSE
         WRITE (6, *) 'ITJPOR: Argument ICSF is out of range.'
         STOP
      ENDIF
!
      RETURN
      END FUNCTION ITJPOR
