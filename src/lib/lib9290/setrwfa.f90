!***********************************************************************
!                                                                      *
      SUBROUTINE SETRWFA(NAME) 
!                                                                      *
!   Open, check, load data from and close the  .rwf  file.             *
!                                                                      *
!   Call(s) to: [LIB92]: LODRWF, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified    by Xinghong He            Last revision: 09 Jul 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:42   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE IOUNIT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      USE lodrwf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = *), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, IOS 
      CHARACTER :: G92RWF*6 


      CALL OPENFL (23, NAME, 'UNFORMATTED', 'OLD', IERR) 
      IF (IERR == 1) THEN 
         WRITE (ISTDE, *) 'Error when opening', NAME(1:LEN_TRIM(NAME)) 
         STOP  
      ENDIF 
!
!   Check the file; if not as expected, stop.
!
      READ (23, IOSTAT=IOS) G92RWF 
      IF (IOS/=0 .OR. G92RWF/='G92RWF') THEN 
         WRITE (ISTDE, *) 'This is not a Radial WaveFunction File;' 
         CLOSE(23) 
         STOP  
      ENDIF 
!
!   Attempt to load the radial wavefunctions; if this fails, stop
!
      CALL LODRWF (IERR) 
 
      IF (IERR /= 0) THEN 
         WRITE (ISTDE, *) 'Radial wavefunctions defined in CSL file', &
            ' not found in Radial WaveFunction File' 
         CLOSE(23) 
         STOP  
      ENDIF 
 
      CLOSE(23) 
 
      RETURN  
      END SUBROUTINE SETRWFA 
