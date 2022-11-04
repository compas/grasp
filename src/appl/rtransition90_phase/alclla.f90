!***********************************************************************
!                                                                      *
      SUBROUTINE ALCLLA(IBEG, ILAB, ILAST, ILEFT, IPTCSF, IRIGHT, LBLINT, &
         LLDIM, IMODE)
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   local to SUBROUTINE TRSORT.                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE memory_man
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LLDIM
      INTEGER , INTENT(IN) :: IMODE
      INTEGER, DIMENSION(:), pointer :: IBEG, ILAB, ILAST, ILEFT, &
                                        IPTCSF, IRIGHT, LBLINT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEWSIZ
!-----------------------------------------------
!
!
      SELECT CASE (IMODE)
      CASE (1)
!
!   Initial array dimension
!
         LLDIM = 64
!
!   Allocate storage for arrays
!
         CALL ALLOC (IBEG, LLDIM, 'IBEG', 'ALCLLA')
         CALL ALLOC (ILAB, LLDIM, 'ILAB', 'ALCLLA')
         CALL ALLOC (ILAST, LLDIM, 'ILAST', 'ALCLLA')
         CALL ALLOC (ILEFT, LLDIM, 'ILEFT', 'ALCLLA')
         CALL ALLOC (IPTCSF, LLDIM, 'IPTCSF', 'ALCLLA')
         CALL ALLOC (IRIGHT, LLDIM, 'IRIGHT', 'ALCLLA')
         CALL ALLOC (LBLINT, LLDIM, 'LBLINT', 'ALCLLA')
!
      CASE (2)
!
!   Double the allocation of storage for the arrays
!
         NEWSIZ = 2*LLDIM
!
         CALL RALLOC (IBEG, NEWSIZ, 'IBEG', 'ALCLLA')
         CALL RALLOC (ILAB, NEWSIZ, 'ILAB', 'ALCLLA')
         CALL RALLOC (ILAST, NEWSIZ, 'ILAST', 'ALCLLA')
         CALL RALLOC (ILEFT, NEWSIZ, 'ILEFT', 'ALCLLA')
         CALL RALLOC (IPTCSF, NEWSIZ, 'IPTCSF', 'ALCLLA')
         CALL RALLOC (IRIGHT, NEWSIZ, 'IRIGHT', 'ALCLLA')
         CALL RALLOC (LBLINT, NEWSIZ, 'LBLINT', 'ALCLLA')
!
         LLDIM = NEWSIZ
!
      CASE (3)
!
!   Deallocate the storage for the arrays
!
         CALL DALLOC (IBEG, 'IBEG', 'ALCLLA')
         CALL DALLOC (ILAB, 'ILAB', 'ALCLLA')
         CALL DALLOC (ILAST, 'ILAST', 'ALCLLA')
         CALL DALLOC (ILEFT, 'ILEFT', 'ALCLLA')
         CALL DALLOC (IPTCSF, 'IPTCSF', 'ALCLLA')
         CALL DALLOC (IRIGHT, 'IRIGHT', 'ALCLLA')
         CALL DALLOC (LBLINT, 'LBLINT', 'ALCLLA')
!
      CASE DEFAULT
!
         WRITE (6, *) 'ALCLLA: Invalid argument IMODE = ', IMODE
         STOP
!
      END SELECT
!
      RETURN
      END SUBROUTINE ALCLLA
