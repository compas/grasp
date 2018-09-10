!***********************************************************************
!                                                                      *
      SUBROUTINE RADGRD 
!                                                                      *
!   This routine sets up the radial grid  R  and the associated arr-   *
!   ays  RP  and  RPOR  in the COMMON block  /GRID/. Different grids   *
!   are generated depending on whether or not  HP  is  zero .          *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:19   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE DEBUG_C 
      USE DEF_C, ONLY: PRECIS  
      USE GRID_C 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP10, I, NB2, NROWS, II, II1, II2 
      REAL(DOUBLE) :: EPH, ETT, ETTM1, EPSLON, A, RLAST, REST, T, RESTS, FOFR, &
         FPRI, DELR 
!-----------------------------------------------
!
!
!   RPOR(1) is never used in the program: it is arbitrarily
!   set to zero
!
      NP10 = N + 10 
      R(1) = 0.0D00 
      RPOR(1) = 0.0D00 
!
!   Now set up the grids
!
      IF (HP == 0.0D00) THEN 
!        default comes here
!
!   Exponential grid if HP is zero
!
!   Initializations
!
         RP(1) = RNT 
         EPH = EXP(H) 
         ETT = 1.0D00 
!
!   Set up the arrays R, RP, RPOR
!
         DO I = 2, NP10 
            ETT = EPH*ETT 
            ETTM1 = ETT - 1.0D00 
            R(I) = RNT*ETTM1 
            RP(I) = RNT*ETT 
            RPOR(I) = ETT/ETTM1 
         END DO 
!
      ELSE 
!
!   Asymptotically-linear exponential grid otherwise:
!
!   Initializations
!
         EPSLON = 1.0D03*PRECIS 
         A = H/HP 
         RP(1) = RNT/(A*RNT + 1.0D00) 
         RLAST = 0.0D00 
         REST = 0.0D00 
!
!   Set up the arrays R, RP, RPOR
!
         DO I = 2, NP10 
!
            T = H*DBLE(I - 1) 
!
!   Solve the implicit equation for R using the Newton-Raphson
!   method
!
    2       CONTINUE 
            RESTS = REST + RNT 
            FOFR = LOG(RESTS/RNT) + A*REST - T 
            FPRI = RESTS/(A*RESTS + 1.0D00) 
            DELR = -FOFR*FPRI 
            REST = RLAST + DELR 
!
            IF (ABS(DELR/REST) < EPSLON) THEN 
               R(I) = REST 
               RESTS = REST + RNT 
               FPRI = RESTS/(A*RESTS + 1.0D00) 
               RP(I) = FPRI 
               RPOR(I) = FPRI/REST 
            ELSE 
               RLAST = REST 
               GO TO 2 
            ENDIF 
!
         END DO 
!
      ENDIF 
!
!   Debug printout
!
      IF (LDBPR(1)) THEN 
         WRITE (99, 300) 
         NB2 = N/2 
         IF (2*NB2 == N) THEN 
            NROWS = NB2 
         ELSE 
            NROWS = NB2 + 1 
         ENDIF 
         DO II = 1, NROWS 
            II1 = II 
            II2 = II1 + NROWS 
            IF (II2 <= N) THEN 
               WRITE (99, 301) R(II1), RP(II1), RPOR(II1), R(II2), RP(II2), &
                  RPOR(II2) 
            ELSE IF (II1 <= N) THEN 
               WRITE (99, 301) R(II1), RP(II1), RPOR(II1) 
            ENDIF 
         END DO 
      ENDIF 
!
      RETURN  
!
  300 FORMAT(/,'From SUBROUTINE RADGRD:'/,2(&
         ' -------- r -------- -------- r'' -------',' ------- r''/r ------')) 
  301 FORMAT(1P,6(1X,1D19.12)) 
      RETURN  
!
      END SUBROUTINE RADGRD 
