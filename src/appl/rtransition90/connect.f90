!***********************************************************************
!                                                                      *
      SUBROUTINE CONNECT 
!                                                                      *
!   The position of an orbital in the merged list is connected to      *
!   the positions in the initial and final state lists                 *
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
      USE vast_kind_param, ONLY:  DOUBLE 
      USE biorb_C
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J 
!-----------------------------------------------
!
!ww      INTEGER PNTRIQ
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     Initialize
!
      NNII(:NW) = 0.D0 
      NNFF(:NW) = 0.D0 
!
!   Loop over the orbitals in the merged list
!
      DO I = 1, NW 
         DO J = 1, NWII 
            IF (NP(I)/=NPII(J) .OR. NAK(I)/=NAKII(J)) CYCLE  
            NNII(I) = J 
         END DO 
         DO J = 1, NWFF 
            IF (NP(I)/=NPFF(J) .OR. NAK(I)/=NAKFF(J)) CYCLE  
            NNFF(I) = J 
         END DO 
      END DO 
 
      RETURN  
      END SUBROUTINE CONNECT 
