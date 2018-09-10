!***********************************************************************
!                                                                      *
      SUBROUTINE BRINT6 (IA,IB,IC,ID,NU,TEGRAL)
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
      IF (.NOT. FIRST(6)) THEN
!
!   This branch is executed on all entries except the first
!
         IF (INDEX .GT. INDTP6(NTPI(6))) THEN
!
!   The index is greater than the largest stored
!
            FOUND = .FALSE.
            LOC = NTPI(6)
!
         ELSEIF (INDEX .LT. INDTP6(1)) THEN
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
            JU = NTPI(6)
            JL = 1
    1       IF (JU-JL .GT. 1) THEN
               JM = (JU+JL)/2
               IF (INDTP6(JM) .GT. INDEX) THEN
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
               IF (INDEX .EQ. INDTP6(JU)) THEN
!
                  FOUND = .TRUE.
                  LOC = JU
!
               ELSEIF (INDEX .EQ. INDTP6(JL)) THEN
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
            TEGRAL = VALTP6(LOC)
!
         ELSE
!
!   Index not found; compute the integral
!
            TEGRAL = BRINTF (6,IA,IB,IC,ID,NU)
!
!   Increment the integral counter
!
            NTPI(6) = NTPI(6)+1
!
!   Increase array length by half the present length if the latter
!   is inadequate to store another pair
!
!
            IF (NTPI(6) .GT. NDTPA(6)) THEN
               NEWSIZ = NDTPA(6)+NDTPA(6)/2
               CALL RALLOC (INDTP6,NEWSIZ,'INDTP6', 'BRINT6')
               CALL RALLOC (VALTP6,NEWSIZ,'VALTP6', 'BRINT6')
               NDTPA(6)= NEWSIZ
            ENDIF
            DO 4 I = NTPI(6),LOC+2,-1
               INDTP6(I) = INDTP6(I-1)
               VALTP6(I) = VALTP6(I-1)
    4       CONTINUE
!
!   Put the new index and value into storage
!
            INDTP6(LOC+1) = INDEX
            VALTP6(LOC+1) = TEGRAL
!
         ENDIF
!
      ELSE
!
!   This branch is executed only once per type of integral
!
         FIRST(6) = .FALSE.
!
!   Designate the initial storage for arrays INDTPx and VALTPx;
!   Array NDTPA stores the array dimensions
!
         NDTPA(6) = 10000
!
!   Compute the integral's value
!
         TEGRAL = BRINTF (6,IA,IB,IC,ID,NU)
!
!   Initialise the integral counter
!
         NTPI(6) = 1
!
!   Store the integral and its value
!
         CALL ALLOC (INDTP6,NDTPA(6),'INDTP6', 'BRINT6')
         CALL ALLOC (VALTP6,NDTPA(6),'VALTP6', 'BRINT6')

         INDTP6(1) = INDEX
         VALTP6(1) = TEGRAL
!
      ENDIF
!
      RETURN
      END
