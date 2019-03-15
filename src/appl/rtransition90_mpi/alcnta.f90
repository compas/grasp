!***********************************************************************
!                                                                      *
      SUBROUTINE ALCNTA(ISLDR, ISLDR1, XSLDR, NTDIM, IMODE)
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   used in OSCL.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
!                                                                      *
!   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
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
      INTEGER  :: NTDIM
      INTEGER , INTENT(IN) :: IMODE
      INTEGER, DIMENSION(:), POINTER :: ISLDR, ISLDR1
      REAL(DOUBLE), DIMENSION(:), POINTER :: XSLDR
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
         NTDIM = 1
!
!   Allocate storage for arrays
!
         CALL ALLOC (ISLDR, NTDIM, 'ISLDR', 'ALCNTA' )
         CALL ALLOC (ISLDR1, NTDIM,'ISLDR1', 'ALCNTA' )
         CALL ALLOC (XSLDR, NTDIM, 'XSLDR', 'ALCNTA' )
!
      CASE (2)
!
!   Double the allocation of storage for the arrays
!
         NEWSIZ = 2*NTDIM
!
         CALL RALLOC (ISLDR, NEWSIZ, 'ISLDR', 'ALCNTA' )
         CALL RALLOC (ISLDR1, NEWSIZ, 'ISLDR1', 'ALCNTA' )
         CALL RALLOC (XSLDR, NEWSIZ, 'XSLDR', 'ALCNTA' )
!
         NTDIM = NEWSIZ
!
      CASE (3)
!
!   Deallocate the storage for the arrays
!
         CALL DALLOC (ISLDR, 'ISLDR', 'ALCNTA' )
         CALL DALLOC (ISLDR1, 'ISLDR1', 'ALCNTA' )
         CALL DALLOC (XSLDR,  'XSLDR', 'ALCNTA' )
!
      CASE DEFAULT
!
         WRITE (6, *) 'ALCNTA: Invalid argument IMODE = ', IMODE
         STOP
!
      END SELECT
!
      RETURN
      END SUBROUTINE ALCNTA
