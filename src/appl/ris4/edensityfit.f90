!***********************************************************************
!                                                                      *
SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NRNUC,F,RHO,RES)
!SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NPARFIT,                  &
!                                                DRMS,P,F,RHO,RES,NRNUC)

!     Fits polynomial b1 + b2r^2 b3r^3 + b4r^4 to (r,rho) electron     *
!     density data using least squares method                          *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE def_C,            ONLY: PI,FMTOAU,CCMS,AUCM
      USE parameter_def,    ONLY: NNNP, NNN1
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNN1) :: XVEC
      REAL(DOUBLE), DIMENSION(NNNP) :: YVEC
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
      REAL(DOUBLE), DIMENSION(5) :: P
      REAL(DOUBLE), DIMENSION(6) :: F
      REAL(DOUBLE), DIMENSION(2) :: PAR
      REAL(DOUBLE) :: Z, DRMS, RHO, RES
      INTEGER      :: NRNUC, NPARFIT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNN1) :: X
      REAL(DOUBLE), DIMENSION(NNNP) :: Y, W, DY
!GG      REAL(DOUBLE), DIMENSION(NNNP) :: X, Y, W, DY
      REAL(DOUBLE), DIMENSION(3,3)  :: M, MI, C, CI
      REAL(DOUBLE), DIMENSION(4)    :: DR, DR2
      REAL(DOUBLE), DIMENSION(3)    :: RM, PM, B, DF2
      REAL(DOUBLE), DIMENSION(2)    :: PARF, PARF2
      REAL(DOUBLE) :: MDET, CDET, AU2FM, NORM, CONST, CONST2, A0
      REAL(DOUBLE) :: A1, A2, FDSUM, FO90, R2,DR4,DR6, DE, DEN, dDEN
      REAL(DOUBLE) :: dDEN2, DX2, DY2, X2, Y2
      INTEGER :: I, INUC, NMIN, NMAX, NR
!-----------------------------------------------
!
      AU2FM = 1.D0/FMTOAU                ! Bohr radius in fm
      CONST = 1.D-9*CCMS*AUCM*AU2FM      ! To get electronic factors in GHz

!     COPY ARRAYS
      X(:) = XVEC(:)
      Y(:) = YVEC(:)

      ! DETERMINE FIRST POINT BEYOND FO90 IN GRID CALLED NR
      ! DETERMINE FIRST POINT WHERE GRID IS RELIABLE CALLED NMIN
      ! SINCE FIRST FEW DATA POINTS IN DENSITY ARE NOT RELIABLE DUE TO DIVISION WITH SMALL
      ! R^2 VALUES (STAGGERING IS SEEN), NMIN IS THE FIRST POINT TO BE USED IN THE SUBSEQUENT LEAST SQUARES FIT

      NMIN = 10

!     SETS NMAX IN FIT TO NR
      NMAX = NRNUC

!     DETERMINE RHO:
!     BELOW PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4 IS FITTED TROUGH DATA POINTS
!     (X(NMIN),Y(NMIN)), (X(NMIN+2),Y(NMIN+2)), (X(NMIN+4),Y(NMIN+4)).
!     WE HAVE: RHO(0) = PM(1)

      RM(1) = Y(NMIN)
      RM(2) = Y(NMIN+2)
      RM(3) = Y(NMIN+4)

      M(1,1) = 1.0d0
      M(1,2) = X(NMIN)**2.0d0
      M(1,3) = X(NMIN)**4.0d0
      M(2,1) = 1.0d0
      M(2,2) = X(NMIN+2)**2.0d0
      M(2,3) = X(NMIN+2)**4.0d0
      M(3,1) = 1.d0
      M(3,2) = X(NMIN+4)**2.0d0
      M(3,3) = X(NMIN+4)**4.0d0
      MDET = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-                &
           M(1,2)*M(2,1)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+                  &
           M(1,3)*M(2,1)*M(3,2)-M(1,3)*M(2,2)*M(3,1)

! Determine inverse of C matrix
      MI(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/MDET
      MI(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/MDET
      MI(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/MDET

      MI(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/MDET
      MI(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/MDET
      MI(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/MDET

      MI(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/MDET
      MI(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/MDET
      MI(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/MDET

      ! Determine parameters
      PM(1) = MI(1,1)*RM(1)+MI(1,2)*RM(2)+MI(1,3)*RM(3)
      PM(2) = MI(2,1)*RM(1)+MI(2,2)*RM(2)+MI(2,3)*RM(3)
      PM(3) = MI(3,1)*RM(1)+MI(3,2)*RM(2)+MI(3,3)*RM(3)

      ! Finallay RHO(0) in au^{-3} is determined
      RHO = PM(1)

      ! START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS
      ! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))

      NORM = -(Y(NMAX)-RHO)/AU2FM**3.0d0

      DO I=NMIN,NMAX
         X(I) = AU2FM*X(I)
         Y(I) = (Y(I)-RHO)/AU2FM**3.0d0
         Y(I) = Y(I)/NORM
         W(I) = X(I)                        ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
      END DO

!     THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
!     P(2)*X(I)**2 + P(3)*X(I)**4 + P(4)*X(I)**6
!     DETERMINE B_l AND C_{kl} MATRIX ELEMENTS
      B(:) = 0.0d0
      C(:,:) = 0.0d0
      DO I=NMIN,NMAX
         B(1) = B(1) + Y(I)*X(I)**2.0d0*W(I)
         B(2) = B(2) + Y(I)*X(I)**4.0d0*W(I)
         B(3) = B(3) + Y(I)*X(I)**6.0d0*W(I)
         C(1,1) = C(1,1) + X(I)**4.0d0*W(I)
         C(2,2) = C(2,2) + X(I)**8.0d0*W(I)
         C(3,3) = C(3,3) + X(I)**12.0d0*W(I)
         C(1,2) = C(1,2) + X(I)**6.0d0*W(I)
         C(1,3) = C(1,3) + X(I)**8.0d0*W(I)
         C(2,3) = C(2,3) + X(I)**10.0d0*W(I)
      END DO
      C(2,1) = C(1,2)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)

!     COMPUTE DETERMINANT
      CDET = C(1,1)*C(2,2)*C(3,3)-C(1,1)*C(2,3)*C(3,2)- &
           C(1,2)*C(2,1)*C(3,3)+C(1,2)*C(2,3)*C(3,1)+ &
           C(1,3)*C(2,1)*C(3,2)-C(1,3)*C(2,2)*C(3,1)

!     COMPUTE INVERSE OF C MATRIX
      CI(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))/CDET
      CI(1,2) = (C(1,3)*C(3,2)-C(1,2)*C(3,3))/CDET
      CI(1,3) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))/CDET

      CI(2,1) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))/CDET
      CI(2,2) = (C(1,1)*C(3,3)-C(1,3)*C(3,1))/CDET
      CI(2,3) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))/CDET

      CI(3,1) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))/CDET
      CI(3,2) = (C(1,2)*C(3,1)-C(1,1)*C(3,2))/CDET
      CI(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))/CDET

!     DETERMINE FITTING PARAMETERS
      P(2) = CI(1,1)*B(1)+CI(1,2)*B(2)+CI(1,3)*B(3)
      P(3) = CI(2,1)*B(1)+CI(2,2)*B(2)+CI(2,3)*B(3)
      P(4) = CI(3,1)*B(1)+CI(3,2)*B(2)+CI(3,3)*B(3)

      P(1) = RHO/AU2FM**3.0d0
      P(2) = P(2)*NORM
      P(3) = P(3)*NORM
      P(4) = P(4)*NORM

!     RESULTING ELECTRONIC FACTORS F(N) IN GHZ fm^{-2N}
      F(1) = 2.0d0*PI/3.0d0*Z*P(1)*CONST
      F(2) = 2.0d0*PI/10.0d0*Z*P(2)*CONST
      F(3) = 2.0d0*PI/21.0d0*Z*P(3)*CONST
      F(4) = 2.0d0*PI/36.0d0*Z*P(4)*CONST

!    DETERMINE AVERAGE POINT DISCREPANCY PARAMETER
      RES = 0.0d0
      DO I=NMIN,NMAX
         Y(I) = Y(I)*AU2FM**3.0d0*NORM+RHO
         RES = RES + (Y(I)-RHO &
              -AU2FM**3.0d0*(P(2)*X(I)**2.0d0+P(3)*X(I)**4.0d0 &
              +P(4)*X(I)**6.0d0))**2.0d0
      END DO

      RES = sqrt(RES/(NMAX-NMIN+1))
      RES = RES/RHO*1000.0d0    ! In per mille of RHO

!     BELOW WE APPROXIMATE THE IS ENERGY. DIFFERENCE IN RADIAL MOMENTS ARE CALUCLATED BETWEEN THE
!     CURRENT DISTRIBTUION (PARF) AND WITH A DISTRIBUTION WITH A R^2 1 FM LARGER. THE DIFFERNCE
!     IN RADIAL MOMENTS ARE THEN MULTIPLIED WITH THE ELECTRONIC FACTORS.

      PARF(1) = PAR(1)*AU2FM    ! 50% FALL OFF RADIUS C
      PARF(2) = PAR(2)*AU2FM    ! SKIN THICKNESS a. 10% TO 90% FALL OFF DISTANCE t = 4.0*ln(3)*a

      !DR2(1) = 1.0d0
      !PARF2(1) =  SQRT(1.d0+PARF(1)**2.d0)
      DR2(2) = 10.d0/7.d0*PARF(1)**2.d0 + 30.d0/7.d0*(PI*PARF(2))**2.d0
      DR2(3) = 5.d0/3.d0*PARF(1)**4.d0 + &
           110.d0/9.d0*(PI*PARF(2)*PARF(1))**2.d0 + &
           239.d0/9.d0*(PI*PARF(2))**4.d0
      DR2(4) = 20.d0/11.d0*PARF(1)**6.d0 + &
           780.d0/33.d0*(PI*PARF(2))**2.d0*PARF(1)**4.d0 + &
           4100.d0/33.d0*(PI*PARF(2))**4.d0*PARF(1)**2.d0 + &
           8180.d0/33.d0*(PI*PARF(2))**6.d0

      FDSUM = F(1)
      DO I=2,4
         FDSUM = FDSUM + F(I)*DR2(I)
      END DO
      F(5) = FDSUM              ! IS ENERGY IN GHZ / FM^2

      DR2(2) = 25.d0/21.d0
      DR2(3) = 25.d0/9.d0*PARF(1)**2.d0 + &
           275.d0/27.d0*(PI*PARF(2))**2.d0
      DR2(4) = 50.d0/11.d0*PARF(1)**4.d0 +  &
           1300.d0/33.d0*(PI*PARF(2)*PARF(1))**2.d0 + &
           10250.d0/99.d0*(PI*PARF(2))**4.d0

      FDSUM = 0.d0
      DO I=2,4
         FDSUM = FDSUM + F(I)*DR2(I)
      END DO
      F(6) = FDSUM              ! IS ENERGY IN GHZ / FM^2

      RETURN

      END SUBROUTINE EDENSITYFIT
