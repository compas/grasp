!***********************************************************************
!                                                                      *
      INTEGER FUNCTION ISPAR (ICSF) 
!                                                                      *
!   ISPAR is the value of P for CSF number ICSF.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 02 Nov 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:41   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE parameter_def, ONLY: NNNW
      USE STAT_C,        ONLY: JCUPA
      USE IOUNIT_C,      ONLY: ISTDE 
      USE ORB_C,         ONLY: NCF
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: ICSF 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
      IF (ICSF>=1 .AND. ICSF<=NCF) THEN 
         ispar = jcupa(NNNW,icsf)
         IF (ISPAR > 127) ISPAR = ISPAR - 256 
         ISPAR = SIGN(1,ISPAR) 
      ELSE 
         WRITE (ISTDE, *) 'ISPAR: Argument ICSF is out of range.' 
         WRITE (ISTDE, *) 'ICSF =',ICSF,' NCF =',NCF 
         STOP  
      ENDIF 
!
      RETURN  
      END FUNCTION ISPAR 
