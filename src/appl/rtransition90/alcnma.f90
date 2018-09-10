!***********************************************************************
!                                                                      *
      SUBROUTINE ALCNMA(IPTR, ISLDR, ISLDR1, XSLDR, NMDIM, IMODE) 
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
      INTEGER  :: NMDIM 
      INTEGER , INTENT(IN) :: IMODE 
      INTEGER, DIMENSION(:), pointer :: IPTR, ISLDR, ISLDR1
      REAL(DOUBLE), DIMENSION(:), pointer :: XSLDR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEWSIZ 
!-----------------------------------------------
!
      SELECT CASE (IMODE)  
      CASE (1)  
!
!   Initial array dimension
!
         NMDIM = 64 
!
!   Allocate storage for arrays
!
         CALL ALLOC (IPTR, NMDIM, 'IPTR', 'ALCNMA') 
         CALL ALLOC (ISLDR, NMDIM, 'ISLDR', 'ALCNMA') 
         CALL ALLOC (ISLDR1, NMDIM, 'ISLDR1', 'ALCNMA') 
         CALL ALLOC (XSLDR, NMDIM, 'XSLDR', 'ALCNMA') 
!
      CASE (2)  
!
!   Double the allocation of storage for the arrays
!
         NEWSIZ = 2*NMDIM 
!
         CALL RALLOC (IPTR, NEWSIZ, 'IPTR', 'ALCNMA') 
         CALL RALLOC (ISLDR, NEWSIZ, 'ISLDR', 'ALCNMA') 
         CALL RALLOC (ISLDR1, NEWSIZ, 'ISLDR1', 'ALCNMA') 
         CALL RALLOC (XSLDR, NEWSIZ, 'XSLDR', 'ALCNMA') 
!
         NMDIM = NEWSIZ 
!
      CASE (3)  
!
!   Deallocate the storage for the arrays
!
         CALL DALLOC (IPTR,  'IPTR', 'ALCNMA') 
         CALL DALLOC (ISLDR, 'ISLDR', 'ALCNMA') 
         CALL DALLOC (ISLDR1, 'ISLDR1', 'ALCNMA') 
         CALL DALLOC (XSLDR,  'XSLDR', 'ALCNMA') 
!
      CASE DEFAULT 
!
         WRITE (6, *) 'ALCNMA: Invalid argument IMODE = ', IMODE 
         STOP  
!
      END SELECT 
!
      RETURN  
      END SUBROUTINE ALCNMA 
