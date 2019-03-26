!***********************************************************************
!                                                                      *
      SUBROUTINE ZKF(K, I, J)
!                                                                      *
!   This subroutine evaluates Hartree Z-functionals:                   *
!                                                                      *
!              (k)                     k                               *
!             Z   [f(r);r] =  I ( (s/r)   f(s) ; 0 - r )               *
!                                                                      *
!   where  I ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
!   range  in  s .    The Z-functional is tabulated in  COMMON/TATB/   *
!   in array  TB . The f-function is assumed tabulted in array  TA .   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNN1
      USE cnc_C, ONLY: cnc5c
      USE grid_C, ONLY: n, r, rp
      USE ncc_C
      USE tatb_C, ZK=>TB
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: K
      INTEGER  :: I
      INTEGER  :: J    !!!! Arument not referenced !!!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: II, MTPP1, MTPP3, MTPP4, KK
      REAL(DOUBLE), DIMENSION(NNN1) :: RHOP, RTTK, TEMP
      REAL(DOUBLE) :: SUM, ZKLIM
!-----------------------------------------------
!
      IF (K > 0) THEN
         RTTK(2:N) = R(2:N)**K
      ENDIF
!
!   MTP is fed in through  COMMON/TATB/
!
      MTPP1 = MTP + 1
      MTPP3 = MTP + 3
      MTPP4 = MTP + 4
!
!   Compute  RP(S)*F(S)  and store it in  RHOP
!
      RHOP(2:MTP) = RP(2:MTP)*TA(2:MTP)
!
!   Fill array TEMP with r**k * RHOP
!
      TEMP(1) = 0.0D00
      IF (K == 0) THEN
         TEMP(2:MTP) = RHOP(2:MTP)
      ELSE
         TEMP(2:MTP) = RTTK(2:MTP)*RHOP(2:MTP)
      ENDIF
!
!   Set an additional four points to zero
!
      TEMP(MTPP1:MTPP4) = 0.0D00
!
!                                     k
!   Compute the first few values of  r  * ZK  using semi-open
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
!                         k
!   Compute remainder of r  * ZK: march out to MTP+3
!
      DO II = 5, MTPP3
         ZK(II) = ZK(II-4) + C1*(TEMP(II-4)+TEMP(II)) + C2*(TEMP(II-3)+TEMP(II-&
            1)) + C3*TEMP(II-2)
      END DO
!                                       k   (k)
!   Determine the asymptotic value of  r * Z
!
!   Compute ZK
!
      ZKLIM = ZK(MTPP3)
!
      IF (K == 0) THEN
!
         ZK(MTPP4:N) = ZKLIM
!
      ELSE
!
         ZK(2:MTPP3) = ZK(2:MTPP3)/RTTK(2:MTPP3)
!
         ZK(MTPP4:N) = ZKLIM/RTTK(MTPP4:N)
!
      ENDIF
!
      RETURN
      END SUBROUTINE ZKF
