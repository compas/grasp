!***********************************************************************
!                                                                      *
      SUBROUTINE QUAD(RESULT)
!                                                                      *
!   The argument result is an approximation  to the integral of F(R)   *
!   from  zero  to  infinity,  where the values of RP(I)*F(R(I)) are   *
!   tabulated in the array  TA(I). The integral in the interval zero   *
!   to R(J) is computed by use of an analytical fit                    *
!                                                                      *
!                                SIGMA                                 *
!                      F(R) = A R                                      *
!                                                                      *
!   A five-point  Closed  Newton-Cotes  formula (cf. F B Hildebrand,   *
!   Introduction to Numerical Analysis, second edition, McGraw-Hill,   *
!   New York, 1974, p 93)  is  used  to  compute the integral in the   *
!   interval  R(J:MTP).  The  contribution  from  the  tail  of  the   *
!   function beyond the last  tabular  point  (MTP) is assumed to be   *
!   negligible. The method uses  MTP+3  tabulation points. Array  TA   *
!   should therefore be dimensioned to at least  N+4 .                 *
!                                                                      *
!   No subroutines called.                                             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:18   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEF_C
      USE GRID_C
      USE NCC_C
      USE TATB_C, ONLY: TA, MTP
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(OUT) :: RESULT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTPM1, I, IP1, LOC
      REAL(DOUBLE) :: TAI, TAIP1, QUOTT, FRIP1, FRI, RATIO, RIP1, RI, SIGMA
!-----------------------------------------------
!
!
!   Find first values that will permit computation of exponent
!
      MTPM1 = MTP - 1
      DO I = 2, MTPM1
!
         TAI = TA(I)
         IF (ABS(TAI) <= 0.D0) CYCLE
!
         IP1 = I + 1
         TAIP1 = TA(IP1)
         QUOTT = TAIP1/TAI
!
         IF (QUOTT <= 0.D0) CYCLE
!
!   Exponent from fit
!
         FRIP1 = TAIP1/RP(IP1)
         FRI = TAI/RP(I)
         RATIO = FRIP1/FRI
         RIP1 = R(IP1)
         RI = R(I)
         SIGMA = LOG(RATIO)/LOG(RIP1/RI)
!
!   Analytical integration and error estimate for interval r(1:i)
!
         FRI = RI*FRI
         RESULT = FRI/(SIGMA + 1.D0)
!
!   Set the tail to zero
!
         TA(MTP+1:3+MTP) = 0.D0
!
!   Newton-Cotes quadature for the remainder
!
         RESULT = RESULT + C1*TAI
         RESULT = RESULT + SUM(C2*(TA(IP1:MTP:4)+TA(IP1+2:MTP+2:4))+C3*TA(IP1+1&
            :MTP+1:4)+C4*TA(IP1+3:MTP+3:4))
         IF (MOD(MTP - I,4) == 0) RESULT = RESULT - C1*TA(MTP)
!
!   Test of result's accuracy; `decomment' to activate
!
!              ESTDER = 10.0D 00*RI*FRI
!              RATIO = ABS (ESTDER/RESULT)
!              IF (RATIO .GT. ACCY) PRINT (*,300) RATIO
!
         GO TO 4
!
      END DO
!
!   No value which will permit computation of exponent
!
      RESULT = 0.D0
!
    4 CONTINUE
      RETURN
!
  300 FORMAT(/,'QUAD: Estimated accuracy is ',1P,D10.3,/,&
         ' Decrease RNT or improve input data conditioning to',' ameliorate.'/)
      RETURN
!
      END SUBROUTINE QUAD
