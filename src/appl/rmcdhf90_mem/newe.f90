!***********************************************************************
!                                                                      *
      SUBROUTINE NEWE(J, SGN, NPRIME, MX, DELEPS, FAIL, INV)
!                                                                      *
!   This  subroutine implements Part 2 of Algorithm 7.1 in  C Froese   *
!   Fischer,  Comput  Phys  Rep, 3 (1986) 273-326. (The present code   *
!   actually  implements  the version used in her program MCHF where   *
!   differences occur.)                                                *
!                                                                      *
!   Call(s) to: [RSCF92]: OUTBND.                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last Update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE int_C
      USE orb_C
      USE scf_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outbnd_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
      INTEGER , INTENT(IN) :: NPRIME
      INTEGER , INTENT(IN) :: MX
      INTEGER , INTENT(OUT) :: INV
      REAL(DOUBLE) , INTENT(IN) :: SGN
      REAL(DOUBLE) , INTENT(INOUT) :: DELEPS
      LOGICAL , INTENT(OUT) :: FAIL
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: P02 = 2.0D-02
      REAL(DOUBLE), PARAMETER :: P05 = 5.0D-02
      REAL(DOUBLE), PARAMETER :: P001 = 1.0D-03
      REAL(DOUBLE), PARAMETER :: P00001 = 1.0D-05
      REAL(DOUBLE), PARAMETER :: D2P5 = 2.5D00
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: EPS, ABDELE, DELEBE, DN, DNPRME, DNPN25, ETRY, DELTA
!-----------------------------------------------
!
!
!   Determine if the iterative process has succeeded
!
      IF (SGN>0.0D00 .AND. MX==0) THEN
         FAIL = .FALSE.
      ELSE
         FAIL = .TRUE.
      ENDIF
!
!   Inversion counter
!
      INV = 0
!
!   If unsuccessful, obtain a new estimate of the eigenvalue
!
      IF (FAIL) THEN
!
!   Define quantities used later
!
         EPS = E(J)
         ABDELE = ABS(DELEPS)
         DELEBE = ABS(DELEPS/EPS)
         DN = DBLE(NP(J))
         DNPRME = DBLE(NPRIME)
         DNPN25 = (DNPRME/DN)**D2P5
!
         IF (ABS(MX)==1 .AND. DELEBE>P02 .OR. MX==0 .AND. DELEBE>=P00001 .AND. &
            ABDELE>=P001) THEN
    1       CONTINUE
            ETRY = EPS + DELEPS
            IF (OUTBND(ETRY) .AND. MX/=0) ETRY = EPS*DNPN25
            IF (OUTBND(ETRY)) ETRY = EPS - DELEPS
            IF (OUTBND(ETRY)) THEN
               DELEPS = 0.5D00*DELEPS
               GO TO 1
            ENDIF
         ELSE IF (MX == 0) THEN
            ETRY = EPS
            INV = 1
            P0 = -P0
            P(:MTP0) = -P(:MTP0)
            Q(:MTP0) = -Q(:MTP0)
            FAIL = .FALSE.
         ELSE IF (SGN < 0.0D00) THEN
            ETRY = EPS*DNPN25
            IF (OUTBND(ETRY)) ETRY = 0.5D00*(EPSMAX + EPS)
            IF (OUTBND(ETRY)) ETRY = 0.5D00*(EPSMIN + EPS)
         ELSE IF (MX < 0) THEN
            DELTA = 1.0D00 - EPS/EPSMAX
            EPSMAX = EPS
            IF (DELTA < P05) THEN
               EMAX = EMAX*DNPN25
            ELSE
               EMAX = EPS*DNPN25
            ENDIF
            ETRY = EMAX
            IF (OUTBND(ETRY)) ETRY = 0.5D00*(EPSMAX + EPSMIN)
         ELSE IF (MX > 0) THEN
            DELTA = 1.0D00 - EPSMIN/EPS
            EPSMIN = EPS
            IF (DELTA < P05) THEN
               EMIN = EMIN*DNPN25
            ELSE
               EMIN = EPS*DNPN25
            ENDIF
            ETRY = EMIN
            IF (OUTBND(ETRY)) ETRY = 0.5D00*(EPSMAX + EPSMIN)
         ENDIF
         E(J) = ETRY
      ENDIF
!
      RETURN
      END SUBROUTINE NEWE
