!***********************************************************************
!                                                                      *
      SUBROUTINE LTAB(IS, NQS, KS, IROWS)
!                                                                      *
!   locates rows of possible parents of active shell states for acc-   *
!   essing  NTAB. It is assumed that empty shells have been elimina-   *
!   ted from consideration by SUBROUTINE RKCO.                         *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:49:59   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE TERMS_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, DIMENSION(4), INTENT(IN) :: IS
      INTEGER, DIMENSION(4), INTENT(INOUT) :: NQS
      INTEGER, DIMENSION(4), INTENT(IN) :: KS
      INTEGER, DIMENSION(4), INTENT(INOUT) :: IROWS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: KQ
      INTEGER :: KQ1, KQ2
      INTEGER  :: I
!-----------------------------------------------
!
!
      IF (IS(1) == IS(2)) NQS(1) = NQS(2) - 1
      IF (IS(3) == IS(4)) NQS(3) = NQS(4) - 1
!
      DO I = 1, 4
!
!   Check that input data are consistent
!
         IF (NQS(I)<=0 .OR. NQS(I)>KS(I)) THEN
            WRITE (*, 300) NQS(I), IS(I), KS(I)
            STOP
         ENDIF
!
         KQ1 = NQS(I) - 1
         KQ2 = KS(I) - KQ1
         KQ(I) = MIN(KQ1,KQ2) + 1
         IF (KQ(I) /= 1) THEN
            IROWS(I) = (KS(I)*(KS(I)-2))/8 + KQ(I)
         ELSE
            IROWS(I) = 1
         ENDIF
!
         IF (IROWS(I) <= NROWS) CYCLE
         WRITE (*, 301)
         STOP
!
      END DO
!
      RETURN
!
  300 FORMAT('LTAB: ',1I3,' Electrons in shell ',1I3,' with 2j+1 = ',1I3)
  301 FORMAT('LTAB: Extend COMMON block TERMS')
      RETURN
!
      END SUBROUTINE LTAB
