!***********************************************************************
!                                                                      *
      SUBROUTINE DPBDT(J) 
!-----------------------------------------------
!                                                                      *
!   This subroutine computes H times the derivative, with respect to   *
!   the internal grid, of the large and small components of the wave   *
!   function with index  J .  These  are tabulated, respectively, in   *
!   arrays  TA  and  TB  in  COMMON  block  /TATB/ .                   *
!                                                                      *
!   A  thirteen-point  Lagrange formaula is used for the calculation   *
!   of derivatives.                                                    *
!                                                                      *
!   Written by Farid F Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:21   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE GRID_C 
      USE LIC13_C 
      USE TATB_C, ONLY: TA, TB
      USE WAVE_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K, IROW, LOC 
      REAL(DOUBLE) :: A1, A2, A3, A4, A5, A6, HDPBDT, HDQBDT, AIK 
!-----------------------------------------------
!

!     Replace the equivalence statement
      A1=A(7,1)
      A2=A(7,2)
      A3=A(7,3)
      A4=A(7,4)
      A5=A(7,5)
      A6=A(7,6)
!
!   Compute derivative in three separate regions
!
!   First, points 1 to 6
!
      DO I = 1,6
         HDPBDT = 0.D0
         HDQBDT = 0.D0
!cjb replace loop by SUM
!        DO 1 K = 1,13
!            AIK = A(I,K)
!            HDPBDT = HDPBDT+AIK*PF(K,J)
!            HDQBDT = HDQBDT+AIK*QF(K,J)
!    1    CONTINUE
         HDPBDT = SUM(A(I   ,:)*PF(:13,J))
         HDQBDT = SUM(A(I   ,:)*QF(:13,J))
!cjb
         TA(I) = HDPBDT
         TB(I) = HDQBDT
      END DO

!   Next, points 7 to N-6
!
!   Special treatment for this region because of the symmetry of
!   the differentiation formula
!
      DO I = 7, N - 6 
         TA(I) = A1*(PF(I - 6,J) - PF(I + 6,J)) + A2*(PF(I - 5,J) - PF(I + 5,J)&
            ) + A3*(PF(I - 4,J) - PF(I + 4,J)) + A4*(PF(I - 3,J) - PF(I + 3,J))&
             + A5*(PF(I - 2,J) - PF(I + 2,J)) + A6*(PF(I - 1,J) - PF(I + 1,J)) 
         TB(I) = A1*(QF(I - 6,J) - QF(I + 6,J)) + A2*(QF(I - 5,J) - QF(I + 5,J)&
            ) + A3*(QF(I - 4,J) - QF(I + 4,J)) + A4*(QF(I - 3,J) - QF(I + 3,J))&
             + A5*(QF(I - 2,J) - QF(I + 2,J)) + A6*(QF(I - 1,J) - QF(I + 1,J)) 
      END DO 
!
!   Last, points N-5 to N
!
      DO I = N - 5, N 
         IROW = I - N + 13 
          HDPBDT = SUM(A(IROW,:)*PF(N-12:N,J))
          HDQBDT = SUM(A(IROW,:)*QF(N-12:N,J))
!- git 1/12/07 code below from the new libraries
!         HDPBDT = 0.D0 
!        HDQBDT = 0.D0 
!        DO K = 1, 13 
!           AIK = A(IROW,K) 
!           LOC = N - 13 + K 
!           HDPBDT = HDPBDT + AIK*PF(LOC,J) 
!           HDQBDT = HDQBDT + AIK*QF(LOC,J) 
!        END DO 
!- git 1/12/07
         TA(I) = HDPBDT 
         TB(I) = HDQBDT 
      END DO 
      RETURN  
      END SUBROUTINE DPBDT 
