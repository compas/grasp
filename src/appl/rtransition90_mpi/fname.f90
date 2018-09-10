!***********************************************************************
!                                                                      *
      SUBROUTINE FNAME(NAME) 
!                                                                      *
!   Determines the name of the initial and final states                *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME(2)*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
!-----------------------------------------------
!
!
!   Obtain the names of the initial and final state files
!
    1 CONTINUE 
      WRITE (6, *) ' Name of the Initial state' 
      READ (*, '(A)') NAME(1) 
      WRITE (6, *) ' Name of the Final state' 
      READ (*, '(A)') NAME(2) 
!
      J = INDEX(NAME(1),' ') 
      IF (J == 1) THEN 
         WRITE (6, *) ' Names may not start with blanks' 
         GO TO 1 
      ENDIF 
!
      J = INDEX(NAME(2),' ') 
      IF (J == 1) THEN 
         WRITE (6, *) ' Names may not start with blanks' 
         GO TO 1 
      ENDIF 
 
      RETURN  
      END SUBROUTINE FNAME 
