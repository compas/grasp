!*******************************************************************
!                                                                  *
      SUBROUTINE GRACAH1(I,J,K,L,M,N,RAC)
!                                                                  *
!   SUBROUTINE TO CALCULATE RACAH COEFFICIENTS.                    *
!   THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE.  *
!   WRITTEN BY N. S. SCOTT                                         *
!   Modified by C. Froese Fischer, March 11, 1988 to use table     *
!                                                         look-up  *
!   and                                                            *
!   G. Gaigalas, September 14, 1997                                *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE, TWO
      USE FACTS_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: J
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(DOUBLE) , INTENT(OUT) :: RAC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1,J2,J3,J4,J5,J6,J7,NUMIN,NUMAX,ICOUNT,KK,KI
!-----------------------------------------------
      J1 = I+J+M
      J2 = K+L+M
      J3 = I+K+N
      J4 = J+L+N
      IF (MOD(J1,2) == 0  .AND.  MOD(J2,2) == 0   .AND. &
          MOD(J3,2) == 0  .AND.  MOD(J4,2) == 0 )  THEN
          J1 = J1/2
          J2 = J2/2
          J3 = J3/2
          J4 = J4/2
          IF (MAX(I,J,M) <= J1 .AND.  MAX(K,L,M) <= J2  .AND. &
              MAX(I,K,N) <= J3 .AND.  MAX(J,L,N) <= J4  )  THEN
              J5 = (I+J+K+L)/2
              J6 = (I+L+M+N)/2
              J7 = (J+K+M+N)/2
              NUMIN = MAX(J1, J2, J3, J4) + 1
              NUMAX = MIN(J5, J6, J7)     + 1
              RAC = ONE
              ICOUNT = 0
              DO KK = NUMIN+1,NUMAX
                KI = NUMAX - ICOUNT
                RAC = ONE - (RAC*(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2)))/ &
                         ((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4))
                ICOUNT = ICOUNT+1
              END DO
              RAC = RAC*EXP(                                          &
                    (GAM(NUMIN+1) - GAM(NUMIN-J1) - GAM(NUMIN-J2) -   &
                     GAM(NUMIN-J3) - GAM(NUMIN-J4) - GAM(J5+2-NUMIN)- &
                     GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN)) +               &
                    (GAM(J1+1-I)+GAM(J1+1-J)+GAM(J1+1-M)-GAM(J1+2) +  &
                     GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)-GAM(J2+2) +  &
                     GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2) +  &
                     GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))/TWO)
              IF (MOD(J5+NUMIN,2) == 0) RAC = -RAC
          ELSE
             RAC = ZERO
          END IF
      ELSE
         RAC = ZERO
      END IF
      RETURN
      END
