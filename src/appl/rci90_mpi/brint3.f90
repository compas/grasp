!***********************************************************************
!                                                                      *
      SUBROUTINE BRINT3 (IA,IB,IC,ID,NU,TEGRAL)
!                                                                      *
!   Returns integrals for the transverse photon interaction.           *
!      Integrals are stored in ordered lists. If the integral cannot   *
!   be read from a list, it is computed by calling BRINTF.             *
!                                                                      *
!   Observ that it is not possible to use more than 100 orbitals when  *
!   the Breit interaction is included. If 100 orbitals are used then   *
!   NU <= 19.                                                          *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!               [RCI92]: BRINTF                                        *
!                                                                      *
!   Written by Farid A Parpia               Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man
      USE bilst_C
      USE orb_C,           ONLY: ncf, nw, iqa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE brintf_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ia, ib, ic, id, nu
      REAL(DOUBLE), INTENT(out) :: tegral
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      LOGICAL :: FOUND
      INTEGER :: key, index, loc, ju, jl, jm, newsiz, i
!-----------------------------------------------
!
      KEY = NW + 1
!
!   Compute the integral label
!
      INDEX = (((NU*KEY+IA)*KEY+IB)*KEY+IC)*KEY+ID
!
      IF (.NOT. FIRST(3)) THEN
!
!   This branch is executed on all entries except the first
!
         IF (INDEX .GT. INDTP3(NTPI(3))) THEN
!
!   The index is greater than the largest stored
!
            FOUND = .FALSE.
            LOC = NTPI(3)
!
         ELSEIF (INDEX .LT. INDTP3(1)) THEN
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
            JU = NTPI(3)
            JL = 1
    1       IF (JU-JL .GT. 1) THEN
               JM = (JU+JL)/2
               IF (INDTP3(JM) .GT. INDEX) THEN
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
               IF (INDEX .EQ. INDTP3(JU)) THEN
!
                  FOUND = .TRUE.
                  LOC = JU
!
               ELSEIF (INDEX .EQ. INDTP3(JL)) THEN
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
            TEGRAL = VALTP3(LOC)
!
         ELSE
!
!   Index not found; compute the integral
!
            TEGRAL = BRINTF (3,IA,IB,IC,ID,NU)
!
!   Increment the integral counter
!
            NTPI(3) = NTPI(3)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            IF (NTPI(3) .GT. NDTPA(3)) THEN
               NEWSIZ = NDTPA(3)+NDTPA(3)/2
               CALL RALLOC (INDTP3,NEWSIZ,'INDTP3', 'BRINT3')
               CALL RALLOC (VALTP3,NEWSIZ,'VALTP3', 'BRINT3')
               NDTPA(3)= NEWSIZ
            ENDIF
            DO 4 I = NTPI(3),LOC+2,-1
               INDTP3(I) = INDTP3(I-1)
               VALTP3(I) = VALTP3(I-1)
    4       CONTINUE
!
!   Put the new index and value into storage
!
            INDTP3(LOC+1) = INDEX
            VALTP3(LOC+1) = TEGRAL
!
         ENDIF
!
      ELSE
!
!   This branch is executed only once per type of integral
!
         FIRST(3) = .FALSE.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         NDTPA(3) = 10000
!
!   Compute the integral's value
!
         TEGRAL = BRINTF (3,IA,IB,IC,ID,NU)
!
!   Initialise the integral counter
!
         NTPI(3) = 1
!
!   Store the integral and its value
!
         CALL ALLOC (INDTP3,NDTPA(3),'INDTP3', 'BRINT3')
         CALL ALLOC (VALTP3,NDTPA(3),'VALTP3', 'BRINT3')
         INDTP3(1) = INDEX
         VALTP3(1) = TEGRAL
!
      ENDIF
!
      RETURN
      END
