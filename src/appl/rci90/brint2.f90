!***********************************************************************
!                                                                      *
      SUBROUTINE BRINT2 (IA,IB,IC,ID,K,TEGRAL)
!                                                                      *
!   Returns integrals for the transverse photon interaction.           *
!                                                                      *
!   Written by Per Jonsson                   Octaober 2014             *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer 
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE bilst_C
      USE orb_C,           ONLY: NW
      USE kkstartbreit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ia, ib, ic, id, k
      REAL(DOUBLE), INTENT(out) :: tegral
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: key, index, loc, ju, jl, jm
!-----------------------------------------------

      KEY = NW + 1

      INDEX = ((IA*KEY+IB)*KEY+IC)*KEY+ID
!
      JL = KSTARTBREIT2(K)
      JU = KSTARTBREIT2(K+1) - 1

      IF (INDEX.LT.INDTP2(JL).OR.INDEX.GT.INDTP2(JU)) THEN
        WRITE(*,*) 'Something wrong in brint2'
        STOP
      ENDIF
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 IF (JU-JL .GT. 1) THEN
        JM = (JU+JL)/2
        IF (INDTP2(JM) .GT. INDEX) THEN
          JU = JM
        ELSE
          JL = JM
        ENDIF
        GOTO 1
      ENDIF
!
!   The range is bracketed to the extent possible
!
      IF (INDEX .EQ. INDTP2(JU)) THEN
        LOC = JU
      ELSEIF (INDEX .EQ. INDTP2(JL)) THEN
        LOC = JL
      ELSE
        WRITE(*,*) K,IA,IB,IC,ID,INDEX
        WRITE(*,*) 'Brint2 Integral not found'
        STOP
      ENDIF
!
!   Return the value of the integral
!   from storage

      TEGRAL = VALTP2(LOC)
!
      RETURN
      END SUBROUTINE BRINT2
