!***********************************************************************
!                                                                      *
      SUBROUTINE VPINT (IA,IB,TEGRAL)
!                                                                      *
!   This  routine computes nuclear vacuum polarization integrals.      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!               [RCI92]: VPINTF.                                       *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
     USE vast_kind_param, ONLY: DOUBLE
     USE parameter_def,   ONLY: KEYORB
     USE memory_man
     USE grid_C
     USE vpilst_C, indoei=>indvpi, valoei=>valvpi, FIRST=>FRSTVP
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     USE vpintf_I
     IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: IA, IB
      REAl(DOUBLE), INTENT(OUT) :: TEGRAL
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FOUND
      INTEGER :: iswap, index, loc, ju, jl, jm, newsiz, i
!------------------------------------------------
!
!   Ensure that the indices are in `canonical' order
!
      IF (IA > IB) THEN
         ISWAP = IB
         IB = IA
         IA = ISWAP
      ENDIF
!
!   Compute the composite (packed) index
!
      INDEX = IA*KEY+IB
!
      IF (.NOT. FIRST) THEN
!
!   This branch is executed on all entries except the first
!
         IF (INDEX > INDOEI(NVPI)) THEN
!
!   The index is greater than the largest stored
!
            FOUND = .FALSE.
            LOC = NVPI
!
         ELSEIF (INDEX < INDOEI(1)) THEN
!
!   The index is less than the smallest stored
!
            FOUND = .FALSE.
            LOC = 0
!
         ELSE
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
            JU = NVPI
            JL = 1
!
    1       IF (JU-JL .GT. 1) THEN
!
               JM = (JU+JL)/2
               IF (INDOEI(JM) .GT. INDEX) THEN
                  JU = JM
               ELSE
                  JL = JM
               ENDIF
               GOTO 1
!
            ELSE
!
!   The range is bracketed to the extent possible
!
               IF (INDEX .EQ. INDOEI(JU)) THEN
!
                  FOUND = .TRUE.
                  LOC = JU
!
               ELSEIF (INDEX .EQ. INDOEI(JL)) THEN
!
                  FOUND = .TRUE.
                  LOC = JL
!
               ELSE
!
                  FOUND = .FALSE.
                  LOC = JL
!
               ENDIF
!
            ENDIF
!
         ENDIF
!
         IF (FOUND) THEN
!
!   Found the index in the list; return the value of the integral
!   from storage
!
            TEGRAL = VALOEI(LOC)
!
         ELSE
!
!   Index not found; compute the integral
!
            TEGRAL = VPINTF (IA,IB)
!
!   Increment the integral counter
!
            NVPI = NVPI+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            IF (NVPI .GT. NDVPA) THEN
               NEWSIZ = NDVPA+NDVPA/2
               CALL RALLOC (INDOEI,NEWSIZ,'INDOEI', 'VPINT')
               CALL RALLOC (VALOEI,NEWSIZ,'VALOEI', 'VPINT')
               NDVPA = NEWSIZ
            ENDIF
!
            DO 2 I = NVPI,LOC+2,-1
               INDOEI(I) = INDOEI(I-1)
               VALOEI(I) = VALOEI(I-1)
    2       CONTINUE
!
!   Put the new index and value into storage
!
            INDOEI(LOC+1) = INDEX
            VALOEI(LOC+1) = TEGRAL
!
         ENDIF
!
      ELSE
!
!   This branch is executed only once
!
         FIRST = .FALSE.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDVPA stores the array dimension
!
         NDVPA = 2
         CALL ALLOC (INDOEI,NDVPA,'INDOEI', 'VPINT')
         CALL ALLOC (VALOEI,NDVPA,'VALOEI', 'VPINT')
!
!   Compute the integral's value
!
         TEGRAL = VPINTF (IA,IB)
!
!   Initialise the integral counter
!
         NVPI = 1
!
!   Store the integral and its value
!
         INDOEI(1) = INDEX
         VALOEI(1) = TEGRAL
!
      ENDIF
!
      RETURN
      END SUBROUTINE VPINT
