!***********************************************************************
!                                                                      *
      SUBROUTINE GETS(S, NSHLI, NSHLF)
!                                                                      *
!   This subroutine calculates the initial and final state             *
!   radial overlap matrix                                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE grid_C
      USE tatb_C
      USE orb_C
      USE biorb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      USE wrtmat_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER       :: NSHLI
      INTEGER       :: NSHLF
      REAL(DOUBLE)  :: S(NSHLI,NSHLF)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: I, J, L
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------
!
      DO I = 1, NSHLI
         DO J = 1, NSHLF
!
!   Determine the maximum tabulation point for the integrand
!
            MTP = MIN(MFII(I),MFFF(J))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
            TA(1) = 0.D0
            DO L = 2, MTP
               TA(L) = (PFII(L,I)*PFFF(L,J) + QFII(L,I)*QFFF(L,J))*RP(L)
            END DO
!
!   Perform the quadrature
!
            CALL QUAD (RESULT)
            S(I,J) = RESULT
         END DO
      END DO
!
!   Print out
!
      WRITE (*, *) '********************'
      WRITE (*, *) ' S matrix from GETS'
      WRITE (*, *) '********************'

      CALL WRTMAT (S, NSHLI, NSHLF, NSHLI, NSHLF)
!
      RETURN
      END SUBROUTINE GETS
