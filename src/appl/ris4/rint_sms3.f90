!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINT_SMS3 (I, J)
!                                                                      *
!   The value of RINT_SMS3 is an approximation to:                     *
!                                                                      *
!                                                                      *
!   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
!   Range.                                                             *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by  G. Gaigalas                                            *
!           and E. Gaudamauskas                                        *
!                        Last             revision:  09 October 2009   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE grid_C
      USE tatb_C
      USE wave_C
      USE def_C
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
      MTP = MIN (MF(I),MF(J))
!
!   Value at first tabulation point is arbitrary
!
      TA(1) = 0.0D00
      DO 1 L = 2,MTP
!HFS        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
        TA(L) =  RP(L)*(QF(L,I)*PF(L,J)                                &
              -         PF(L,I)*QF(L,J))/R(L)
    1 CONTINUE
!
!   Perform integration
!
      CALL QUAD (RESULT)
      RINT_SMS3 = RESULT*Z/C
!
      RETURN
      END FUNCTION RINT_SMS3
