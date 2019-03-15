!***********************************************************************
!                                                                      *
      INTEGER FUNCTION JCUPR (LOC, ICSF)
!                                                                      *
!   JCUPR is the 2J+1 value of the LOCth nontrivial intermediate ang-  *
!   ular momentum in CSF  ICSF.                                        *
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
      use orb_C,         ONLY: NWR, NCFR
      use stat_C,        ONLY: JCUPAR
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LOC
      INTEGER  :: ICSF
!-----------------------------------------------
!
      IF (LOC>=1 .AND. LOC<=NWR-1) THEN
         IF (ICSF>=1 .AND. ICSF<=NCFR) THEN
            jcupr = jcupar(loc,icsf)
         ELSE
            WRITE (6, *) 'JCUPR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         WRITE (6, *) 'JCUPR: Argument LOC is out of range.'
         STOP
      ENDIF
!
      RETURN
      END FUNCTION JCUPR
