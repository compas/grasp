!***********************************************************************
!                                                                      *
      SUBROUTINE PACK(IUNPKD, ISUBSH, IPACKD)
!                                                                      *
!   Subshell occupation numbers and all angular momenta 2J+1 are not   *
!   likely to exceed 127 in any application of the GRASP92 suite. It   *
!   is, therefore, inefficient to allocate an entire INTEGER storage   *
!   cell to  any of these quantities --- a single  byte is adequate.   *
!   Up to eight integers of magnitude less than or equal to  127 may   *
!   be stored in one  64-bit cell, four in a 32-bit cell.  This idea   *
!   is implemented in  the present subprogram.  IPACKD is assumed to   *
!   be an INTEGER vector of at least NINT (NW/8) elements for 64-bit   *
!   architectures, and NINT (NW/4) elements for 32-bit architectures.  *
!   ISUBSH  is the subshell sequence number, IUNPKD its value.  LOC1   *
!   is the element number;  LOC2 is one less than the byte number in   *
!   this element.                                                      *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:08   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  BYTE
      USE IOUNIT_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: IUNPKD
      INTEGER, INTENT(IN) :: ISUBSH
      INTEGER(BYTE), DIMENSION(*), INTENT(INOUT) :: IPACKD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
      IF (ABS(IUNPKD) > 127) THEN
         WRITE (ISTDE, *) 'PACK: Argument IUNPKD out of range.'
         STOP
      ENDIF
!
      IF (ISUBSH <= 0) THEN
         WRITE (ISTDE, *) 'PACK: ISUBSH=', ISUBSH, ' less than 1'
         STOP
      ENDIF
!
      IPACKD(ISUBSH) = IUNPKD
!
      RETURN
      END SUBROUTINE PACK
