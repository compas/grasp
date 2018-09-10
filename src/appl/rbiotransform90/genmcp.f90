!***********************************************************************
!                                                                      *
      SUBROUTINE GENMCP(NAME, IC, NTESTG, INPCI) 
!                                                                      *
!   Entry routine for GENMCP. Controls the computation of the          *
!   one-particle coupling coefficients.                                *
!                                                                      *
!   Written by Per Jonsson                                 June 1996   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcslb_I 
      USE factt_I 
      USE mcpout_I 
      USE mcpin_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IC 
      INTEGER  :: NTESTG 
      INTEGER  :: INPCI 
      CHARACTER  :: NAME*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCORE_NOT_USED 
!-----------------------------------------------
!
      WRITE (6, *) 
      WRITE (6, *) 'GENMCP: Execution begins for ', NAME 
!
!   Open, check, load data from, and close the  csl  file
!
      CALL SETCSLB (NAME, NCORE_NOT_USED,3) 
!
!   Set up the table of logarithms of factorials for use by
!   angular modules
!
      CALL FACTT 
!
!   Proceed with the generation of MCP coefficients
!
      CALL MCPOUT (NAME, IC, NTESTG, INPCI) 
      CALL MCPIN (NAME, IC, NTESTG, INPCI) 
!
!   Print completion message
!
      WRITE (6, *) 
!
      RETURN  
      END SUBROUTINE GENMCP 
