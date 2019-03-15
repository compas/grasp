!***********************************************************************
!                                                                      *
      SUBROUTINE ALCNSA(JJA, JJB, HB1, HB2, HC1, HC2, HM1, &
          HM2, LAB, NPTR, NSDIM, IMODE)
!                                                                      *
!   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
!   and deallocates (IMODE = 3) storage for  certain arrays that are   *
!   used in 12_oscl.                                                   *
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
      INTEGER  :: NSDIM
      INTEGER, INTENT(IN) :: IMODE
      INTEGER, DIMENSION(:), pointer :: JJA, JJB, LAB, NPTR
      REAL(DOUBLE), DIMENSION(:), POINTER :: HB1, HB2, HC1, HC2, HM1, HM2
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
         NSDIM = 1
!
!   Allocate storage for arrays
!
         CALL ALLOC (JJA, NSDIM, 'JJA', 'ALCNSA')
         CALL ALLOC (JJB, NSDIM, 'JJB', 'ALCNSA' )
         CALL ALLOC (HB1, NSDIM, 'HB1', 'ALCNSA')
         CALL ALLOC (HB2, NSDIM, 'HB2', 'ALCNSA')
         CALL ALLOC (HC1, NSDIM, 'HC1', 'ALCNSA')
         CALL ALLOC (HC2, NSDIM, 'HC2', 'ALCNSA')
         CALL ALLOC (HM1, NSDIM, 'HM1', 'ALCNSA')
         CALL ALLOC (HM2, NSDIM, 'HM2', 'ALCNSA')
         CALL ALLOC (LAB, NSDIM, 'LAB', 'ALCNSA' )
         CALL ALLOC (NPTR, NSDIM, 'NPTR', 'ALCNSA' )
!
      CASE (2)
!
!   Double the allocation of storage for the arrays
!
         NEWSIZ = 2*NSDIM
!
         CALL RALLOC (JJA, NEWSIZ, 'JJA', 'ALCNSA')
         CALL RALLOC (JJB, NEWSIZ, 'JJB', 'ALCNSA' )
         CALL RALLOC (HB1, NEWSIZ, 'HB1', 'ALCNSA')
         CALL RALLOC (HB2, NEWSIZ, 'HB2', 'ALCNSA')
         CALL RALLOC (HC1, NEWSIZ, 'HC1', 'ALCNSA')
         CALL RALLOC (HC2, NEWSIZ, 'HC2', 'ALCNSA')
         CALL RALLOC (HM1, NEWSIZ, 'HM1', 'ALCNSA')
         CALL RALLOC (HM2, NEWSIZ, 'HM2', 'ALCNSA')
         CALL RALLOC (LAB, NEWSIZ, 'LAB', 'ALCNSA' )
         CALL RALLOC (NPTR, NEWSIZ, 'NPTR', 'ALCNSA' )
!
         NSDIM = NEWSIZ
!
      CASE (3)
!
!   Deallocate the storage for the arrays
!
         CALL DALLOC (JJA, 'JJA', 'ALCNSA')
         CALL DALLOC (JJB, 'JJB', 'ALCNSA' )
         CALL DALLOC (HB1, 'HB1', 'ALCNSA')
         CALL DALLOC (HB2, 'HB2', 'ALCNSA')
         CALL DALLOC (HC1, 'HC1', 'ALCNSA')
         CALL DALLOC (HC2, 'HC2', 'ALCNSA')
         CALL DALLOC (HM1, 'HM1', 'ALCNSA')
         CALL DALLOC (HM2, 'HM2', 'ALCNSA')
         CALL DALLOC (LAB, 'LAB', 'ALCNSA' )
         CALL DALLOC (NPTR, 'NPTR', 'ALCNSA' )
!
      CASE DEFAULT
!
         WRITE (6, *) 'ALCNSA: Invalid argument IMODE = ', IMODE
         STOP
!
      END SELECT
!
      RETURN
      END SUBROUTINE ALCNSA
