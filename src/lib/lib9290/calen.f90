!***********************************************************************
!                                                                      *
      SUBROUTINE CALEN(JTIME, JDATE) 
!                                                                      *
!   Loads the character strings JTIME and JDATE with the time of day   *
!   and the date when called.                                          *
!                                                                      *
!                                         Last revision: 25 Sep 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:38   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: JTIME*10 
      CHARACTER  :: JDATE*8 

      CALL DATE_AND_TIME(JDATE, JTIME)
      
      RETURN  
      END SUBROUTINE CALEN 
