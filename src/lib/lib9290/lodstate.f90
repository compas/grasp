!***********************************************************************
!
      SUBROUTINE LODSTATE(IDBLK) 
!
!   Print block info and ask ASF serial numbers for each block
!
!   Input:
!      nblock, idblk(), ncfblk()
!
!   Output:
!      nevblk(), ncmaxblk()
!      ncmin, iccmin(1:ncmin) -- via items, in common block
!
!    Memories allocated here but not de-allocated are
!        pccmin (via items, deallocated here only if ncmin=0)
!
!  Written by Xinghong He                                  Jul 17 1997
!  Updated by Xinghong He                                  Jun 10 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:30:40   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE memory_man
      USE DEF_C 
      USE hblock_C
      use iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE items_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = 8), DIMENSION(*) :: IDBLK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, JBLOCK, NCF, NCMINOLD, IERR, NTMP 
      CHARACTER :: STR*256 
!-----------------------------------------------
!
 
      WRITE (ISTDE, *) 'There are ', NBLOCK, ' blocks  ', &
         '(block   J/Parity   NCF):' 
      WRITE (ISTDE, '( 4(I3, 1X, A5, I8, 5X) )') (J,IDBLK(J)(1:5),NCFBLK(J),J=1&
         ,NBLOCK) 

      NEVBLK = 0 
      NCMAXBLK = 0 
 
      WRITE (ISTDE, *) 
      WRITE (ISTDE, *) 'Enter ASF serial numbers for each block' 
 
      NCMIN = 0 
  123 CONTINUE 
      DO JBLOCK = 1, NBLOCK 
         NCF = NCFBLK(JBLOCK) 
  234    CONTINUE 
         WRITE (ISTDE, *) 'Block ', JBLOCK, '   ncf = ', NCF, ' id = ', IDBLK(&
            JBLOCK)(1:5) 
 
!        ...Read and parse the list of levels
 
         READ (*, '(A)') STR 
         WRITE(734,'(a)') trim(str) ! write to rscf.log file see, rscf
 
!        ...ICCMIN is allocated and accumulated in items
!        ...ncmin is both input and output parameters to items
         NCMINOLD = NCMIN 
         CALL ITEMS (NCMIN, NCF, STR, IERR) 
         IF (NCMIN == 0) CALL DALLOC (ICCMIN, 'ICCMIN', 'LODSTATE' ) 
         IF (IERR < 0) GO TO 234 
         NEVBLK(JBLOCK) = NCMIN - NCMINOLD 
 
!        ...Determine ncmaxblk
         NTMP = 0 
         DO I = NCMINOLD + 1, NCMIN 
            NTMP = MAX(NTMP,ICCMIN(I)) 
         END DO 
         NCMAXBLK(JBLOCK) = NTMP 
      END DO 
 
      IF (NCMIN == 0) THEN 
         WRITE (ISTDE, *) 'At least one state should be selected' 
         GO TO 123 
      ENDIF 
 
      RETURN  
      END SUBROUTINE LODSTATE 
