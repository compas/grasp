!***********************************************************************
!                                                                      *
      SUBROUTINE IABINT (IA,IB,TEGRAL)
!                                                                      *
!   This routine returns I(ab) integrals.                              *
!      A list of previously computed integrals is maintained. If the   *
!   required  integral is not already  available it is computed  and   *
!   stored in the list.                                                *
!      Note that indices IA and IB may be interchanged on output.      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC, RINTI                          *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE memory_man
      USE coeils_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rinti_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: ia, ib
      REAL(DOUBLE), INTENT(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  FOUND
      INTEGER ::  i, index, loc, ju, jl, jm, key, newsiz, iswap
!-----------------------------------------------
!
!   Ensure that the indices are in `canonical' order
!
      IF (IA .GT. IB) THEN
         ISWAP = IB
         IB = IA
         IA = ISWAP
      ENDIF
!
!   Compute the composite (packed) index
!
      INDEX = IA*KEYORB+IB
!
      IF (.NOT. FRSTCO) THEN
!
!   This branch is executed on all entries except the first
!
         IF (INDEX .GT. INDOEI(NCOEI)) THEN
!
!   The index is greater than the largest stored
!
            FOUND = .FALSE.
            LOC = NCOEI
!
         ELSEIF (INDEX .LT. INDOEI(1)) THEN
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
            JU = NCOEI
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
            TEGRAL = RINTI (IA,IB,0)
!
!   Increment the integral counter
!
            NCOEI = NCOEI+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
            IF (NCOEI .GT. NDCOEA) THEN
               NEWSIZ = NDCOEA+NDCOEA/2
               CALL RALLOC (INDOEI,NEWSIZ,'INDOEI','IABINT' )
               CALL RALLOC (VALOEI,NEWSIZ,'VALOEI','IABINT')
               NDCOEA = NEWSIZ
            ENDIF
!
            DO 2 I = NCOEI,LOC+2,-1
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
         FRSTCO = .FALSE.
!
!   Allocate the initial storage for arrays INDOEI and VALOEI;
!   NDCOEA stores the array dimension
!
         NDCOEA = 2
         CALL ALLOC (INDOEI,NDCOEA,'INDOEI','IABINT' )
         CALL ALLOC (VALOEI,NDCOEA,'VALOEI','IABINT')
!
!   Compute the integral's value
!
         TEGRAL = RINTI (IA,IB,0)
!
!   Initialise the integral counter
!
         NCOEI = 1
!
!   Store the integral and its value
!
         INDOEI(1) = INDEX
         VALOEI(1) = TEGRAL
!
      ENDIF
!
      RETURN
      END SUBROUTINE IABINT
