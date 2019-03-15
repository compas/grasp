!***********************************************************************
!                                                                      *
      SUBROUTINE BRKT
!                                                                      *
!   This subroutine calculates the initial and final state             *
!   radial overlap matrix                                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:56:42   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE grid_C
      USE tatb_C
      USE biorb_C
      USE wave_C, ONLY: mfii, pfii, qfii, mfff, pfff, qfff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, L
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: BRAKET
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------

!
      DO I = 1, NWII
         DO J = 1, NWFF
            IF (NAKII(I) /= NAKFF(J)) CYCLE
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

            BRAKET(I,J) = RESULT

            WRITE (*,9) '<',NPII(I),NHII(I),'|',NPFF(J),NHFF(J),'> =', &
               BRAKET(I,J)
         END DO
      END DO
!
    9 FORMAT(A,I2,A,A,I2,A,A,E20.13)

      RETURN
      END SUBROUTINE BRKT
