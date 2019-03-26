!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINTISO (I, J)
!                                                                      *
!   The value of RINTISO is an approximation to:                       *
!                                                                      *
!                                                                      *
!         I (dv(r)*( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)        *
!                     I     J       I     J                            *
!                                                                      *
!   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
!   Range and dv(r) is the difference between the potentials of two    *
!   Fermi nuclear distributions                                        *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by Per Jonsson             Last revision: 24 Dec 1992      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE grid_C
      USE dvpot_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: L
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------
!
!   Tabulate integrand as required for SUBROUTINE QUAD
!
      MTP = MIN(MF(I),MF(J))
!
!   Value at first tabulation point is arbitrary
!
      TA(1) = 0.0D00
      DO L = 2, MTP
         TA(L) = DV(L)*(PF(L,I)*PF(L,J) + QF(L,I)*QF(L,J))*RP(L)
      END DO
!
!   Perform integration
!
      CALL QUAD (RESULT)
      RINTISO = RESULT
!
      RETURN
      END FUNCTION RINTISO
