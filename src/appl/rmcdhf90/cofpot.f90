      SUBROUTINE COFPOT(EOL, J, NPTS) 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE pote_C
      USE mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setcof_I 
      USE ypot_I 
      USE xpot_I 
      USE lagcon_I 
      USE dacon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
      INTEGER  :: NPTS 
      LOGICAL  :: EOL 
!-----------------------------------------------------------------------

      CALL SETCOF (EOL, J) 
      CALL YPOT (J) 
      CALL XPOT (J) 
      CALL LAGCON (J, NPROCS) 
      CALL DACON 
 
      RETURN  
      END SUBROUTINE COFPOT 
