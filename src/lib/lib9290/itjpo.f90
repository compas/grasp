!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ITJPO (ICSF)
!                                                                      *
!   ITJPO is the value of 2J+1 for CSF number ICSF.                    *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:45   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE STAT_C,        ONLY: JCUPA
      USE IOUNIT_C,      ONLY: ISTDE
      USE orb_C,         ONLY: NCF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: ICSF
!-----------------------------------------------
      IF (ICSF>=1 .AND. ICSF<=NCF) THEN
        itjpo = jcupa(NNNW,icsf)
        IF (ITJPO > 127) ITJPO = 256 - ITJPO
        ITJPO = IABS (ITJPO)
      ELSE
         WRITE (ISTDE, *) 'ITJPO: Argument ICSF is out of range.'
         STOP
      ENDIF
!
      RETURN
      END FUNCTION ITJPO
