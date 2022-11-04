!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ISPARR (ICSF)
!                                                                      *
!   ISPARR is the value of P for CSF number ICSF.                      *
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
      USE ORB_C,         ONLY: NCFR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ICSF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
      IF (ICSF>=1 .AND. ICSF<=NCFR) THEN
         isparr = jcupar(NNNW,icsf)

         IF (ISPARR > 127) ISPARR = ISPARR - 256
         ISPARR = SIGN(1,ISPARR)
      ELSE
         WRITE (6, *) 'ISPARR: Argument ICSF is out of range.'
         STOP
      ENDIF
!
      RETURN
      END FUNCTION ISPARR
