!***********************************************************************
!                                                                      *
      INTEGER FUNCTION INDTPI (ITYPE, I)
!                                                                      *
!   Written by Farid A Parpia             Last revision: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      USE bilst_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ITYPE
      INTEGER  :: I
!-----------------------------------------------
!
      SELECT CASE (ITYPE)
      CASE (1)
         INDTPI = INDTP1(I)
      CASE (2)
         INDTPI = INDTP2(I)
      CASE (3)
         INDTPI = INDTP3(I)
      CASE (4)
         INDTPI = INDTP4(I)
      CASE (5)
         INDTPI = INDTP5(I)
      CASE (6)
         INDTPI = INDTP6(I)
      END SELECT
!
      RETURN
      END FUNCTION INDTPI
