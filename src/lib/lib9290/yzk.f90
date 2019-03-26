!***********************************************************************
!                                                                      *
      SUBROUTINE YZK(K, I, J)
!                                                                      *
!   This subroutine evaluates Hartree Y- and Z-functions:              *
!                                                                      *
!               (K)            (K)           (K)                       *
!              Y   (I,J;r) =  Z   (I,J;r) + W   (I,J;r)                *
!                                                                      *
!   where                                                              *
!                                                                      *
!    (K)                                                               *
!   Z   (I,J;r) =  I ( (s/r)   (P (s)*P (s) + Q (s)*Q (s)) ; 0 - r )   *
!                                I     J       I     J                 *
!                                                                      *
!   and                                                                *
!                                                                      *
!    (K)                    K+1                                        *
!   W   (I,J;r) =  I ( (r/s)   (P (s)*P (s) + Q (s)*Q (s)) ; r -       *
!                                I     J       I     J    INFINITY )   *
!                                                                      *
!   where  I ( G(r,s) ; range )  denotes the integral of G(r,s) over   *
!   range  in  s .  The Y-function is tabulated in  COMMON/TATB/  in   *
!   array  TB , the Z-function in array tA .                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!   Modified by Anders Ynnerman, at Vanderbilt         : 03 Feb 1994   *
!   Modified by Jacek Bieron     at Vanderbilt         : 08 Feb 1994   *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:01   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNN1
      USE CNC_C,           ONLY: CNC5C
      USE DEF_C,           ONLY: ACCY
      USE GRID_C
      USE NCC_C
      USE ORB_C
      USE TATB_C,          ONLY: ZK=>TA, YK=>TB, MTP
      USE WAVE_C,          ONLY: MF , PF, QF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: K
      INTEGER  :: I
      INTEGER  :: J
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KLIMIT = 20
      INTEGER, PARAMETER :: KLIMIT1 = KLIMIT + 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: II, MTPP1, MTPP3, MTPP4, KK, NP4
      REAL(DOUBLE), DIMENSION(NNN1) :: RHOP
      REAL(DOUBLE), DIMENSION(NNN1,KLIMIT) :: RTTK, RTTKM, RTTK1, RTTKM1
      REAL(DOUBLE), DIMENSION(NNN1) :: RM, WK, TEMP
      REAL(DOUBLE) :: SUM, RTMP, ZKLIM, DIF
      LOGICAL, DIMENSION(0:KLIMIT) :: KCALC

      SAVE KCALC, RTTK, RTTKM, RTTK1, RTTKM1, RM
!-----------------------------------------------
!
!      P_OINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
!
!
      DATA KCALC/ KLIMIT1*.FALSE./
      !EQUIVALENCE (TA(1), ZK(1)), (TB(1), YK(1))
      IF (.NOT.KCALC(K)) THEN
         KCALC(K) = .TRUE.
         IF (K > KLIMIT) THEN
            WRITE (6, *) ' increase klimit in yzk.f '
            STOP
         ENDIF
         IF (K > 0) THEN
            RTTK(2:N,K) = R(2:N)**K
            RTTKM(2:N,K) = 1.D0/RTTK(2:N,K)
            RTTK1(2:N,K) = RTTK(2:N,K)*R(2:N)
            RTTKM1(2:N,K) = 1.D0/RTTK1(2:N,K)
            RM(2:N) = 1.D0/R(2:N)
         ELSE
            RM(2:N) = 1.D0/R(2:N)
         ENDIF
      ENDIF
!
!   Determine maximum tabulation point as location beyond which
!   RHOP  (see comment statements below) would be zero; determine
!   other important locations
!
      MTP = MIN(MF(I),MF(J))
      MTPP1 = MTP + 1
      MTPP3 = MTP + 3
      MTPP4 = MTP + 4
!
!   Compute RP(s)*(P (s)*P (s)+Q (s)*Q (s)) and store in RHOP
!                   I     J     I     J
!
      DO II = 2, MTP
         RHOP(II) = RP(II)*(PF(II,I)*PF(II,J) + QF(II,I)*QF(II,J))
      END DO
!
!   Fill array TEMP with R**K * RHOP
!
      TEMP(1) = 0.0D00
      IF (K == 0) THEN
         TEMP(2:MTP) = RHOP(2:MTP)
      ELSE
         TEMP(2:MTP) = RTTK(2:MTP,K)*RHOP(2:MTP)
      ENDIF
!
!   Set an additional four points to zero
!
      TEMP(MTPP1:MTPP4) = 0.0D00
!
!                                     K
!   Compute the first few values of  R  * ZK  using semi-open
!   Newton-Cotes formulae
!
      ZK(1) = 0.0D00
      DO II = 2, 4
         SUM = 0.0D00
         DO KK = 2, 5
            SUM = SUM + CNC5C(KK,II)*TEMP(KK)
         END DO
         ZK(II) = SUM
      END DO
!                         K
!   Compute remainder of R  * ZK: march out to MTP+3; use closed
!   Newton-Cotes formula
!
      DO II = 5, MTPP3
         RTMP = C1*(TEMP(II-4)+TEMP(II)) + C2*(TEMP(II-3)+TEMP(II-1)) + C3*TEMP&
            (II-2)
         ZK(II) = ZK(II-4) + RTMP
      END DO
!
!                                       K   (K)
!   Determine the asymptotic value of  R * Z
!
!                   (0)
!   Correction to  z   : in the manner of  C Froese Fischer,
!   The Hartree-Fock Method for Atoms, John Wiley & Sons,
!   New York, 1977, p 235.
!
      IF (K == 0) THEN
!
         IF (I == J) THEN
            ZKLIM = 1.0D00
         ELSE
            ZKLIM = 0.0D00
         ENDIF
!
         DO KK = MTPP3, MTP, -1
            DIF = ZK(KK) - ZKLIM
            IF (ABS(DIF) <= ACCY) CYCLE
            ZK(KK:KK-((-2-KK)/(-4)-1)*4:(-4)) = ZK(KK:KK-((-2-KK)/(-4)-1)*4:&
               (-4)) - DIF
         END DO
!
      ELSE
!
         ZKLIM = ZK(MTPP3)
!
      ENDIF
!
!   Tabulate  ZK  for entire internal grid
!
      IF (K == 0) THEN
!
         ZK(MTPP4:N) = ZKLIM
!
      ELSE
!
         ZK(2:MTPP3) = ZK(2:MTPP3)*RTTKM(2:MTPP3,K)
!
         ZK(MTPP4:N) = ZKLIM*RTTKM(MTPP4:N,K)
!
      ENDIF
!
!   Start array WK / R**(K+1)
!
      NP4 = N + 4
      WK(NP4:MTPP1:(-1)) = 0.0D00
!
!   Fill array TEMP with RHOP / R**(K+1) ; set TEMP(1) = 0
!   to avoid 0/0 case
!
      TEMP(1) = 0.0D00
      IF (K == 0) THEN
         TEMP(2:MTP) = RHOP(2:MTP)*RM(2:MTP)
      ELSE
         TEMP(2:MTP) = RHOP(2:MTP)*RTTKM1(2:MTP,K)
      ENDIF
!
!   Compute remainder of WK / R**(K+1): march in to the origin
!
      DO II = MTP, 2, -1
         WK(II) = WK(II+4) + C1*(TEMP(II)+TEMP(II+4)) + C2*(TEMP(II+1)+TEMP(II+&
            3)) + C3*TEMP(II+2)
      END DO
      WK(1) = 0.0D00
!
!   Compute WK
!
      IF (K == 0) THEN
         WK(2:MTP) = WK(2:MTP)*R(2:MTP)
      ELSE
         WK(2:MTP) = WK(2:MTP)*RTTK1(2:MTP,K)
      ENDIF
!
!   Assemble solution
!
      YK(1) = 0.0D00
      YK(2:N) = ZK(2:N) + WK(2:N)
!
      RETURN
      END SUBROUTINE YZK
