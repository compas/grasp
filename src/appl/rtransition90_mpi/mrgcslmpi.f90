!***********************************************************************
!                                                                      *
      SUBROUTINE MRGCSL(NAME)
!                                                                      *
!   Entry routine for merging two csl lists                            *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE debug_C
      USE def_C, ONLY:  EMN, IONCTY, NELEC, Z
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ldcsl1_I
      USE ldcsl2_I
      USE merg12_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME(2)*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCORER, NCORE
!-----------------------------------------------
!

!      WRITE (6, *)
!      WRITE (6, *) 'MRGCSL: Execution begins ...'
!
!   Load the first  .csl  file
!
      CALL LDCSL1 (NCORER, NAME(1))
!
!   Load the second  .csl  file
!
      CALL LDCSL2 (NCORE, NAME(2))
!
!   Merge the two  .csl  lists, observe that there may be doublets
!   among the CSF's
!
      CALL MERG12 (NAME, NCORER, NCORE)

      RETURN
      END SUBROUTINE MRGCSL
