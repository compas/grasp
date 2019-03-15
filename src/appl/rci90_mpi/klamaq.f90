!***********************************************************************
!                                                                      *
      SUBROUTINE KLAMAQ(N, KAPPA, Z, FZALFA)
!                                                                      *
!   The function  F (Z*\alpha)  is estimated here. We use the series   *
!   expansion given by  Eqs (1) and (2) and the table of Bethe loga-   *
!   rithms in  S Klarsfeld and A Maquet, Physics Letters  43B (1973)   *
!   201. The vacuum-polarization contribution in Eq (2) is omitted.    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE def_C,           ONLY: C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: KAPPA
      REAL(DOUBLE) , INTENT(IN) :: Z
      REAL(DOUBLE) , INTENT(OUT) :: FZALFA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L, LOC
      REAL(DOUBLE), DIMENSION(36) :: BETHE
      REAL(DOUBLE) :: C401, C402, OVLFAC, BETHEL, TERM, ZALFA, FACTOR
      LOGICAL :: FIRST
!-----------------------------------------------
!
      DATA BETHE/ 2.9841285D00, 2.8117699D00, -0.0300167D00, 2.7676636D00, &
         -0.0381902D00, -0.0052321D00,  2.7498118D00, -0.0419549D00, &
         -0.0067409D00, -0.0017337D00,  2.7408237D00, -0.0440347D00, &
         -0.0076008D00, -0.0022022D00, -0.0007721D00,  2.7356642D00, &
         -0.0453122D00, -0.0081472D00, -0.0025022D00, -0.0009628D00, &
         -0.0004079D00,  2.7324291D00, -0.0461552D00, -0.0085192D00, &
         -0.0027091D00, -0.0010945D00, -0.0004997D00, -0.0002409D00, &
          2.7302673D00, -0.0467413D00, -0.0087850D00, -0.0028591D00, &
         -0.0011904D00, -0.0005665D00, -0.0002904D00, -0.0001539D00/
!
!----------------------------------------------------------------------*
!
      DATA FIRST/ .TRUE./
!
      DATA C401/ 0.0D00/
      DATA C402/ 0.0D00/
      DATA OVLFAC/ 0.0D00/
!
!   Set up the constants
!
      IF (FIRST) THEN
!
         C401 = 11.0D00/24.0D00
         C402 = 3.0D00/8.0D00
         OVLFAC = 4.0D00/3.0D00
!
         FIRST = .FALSE.
!
      ENDIF
!
!   Ensure that the principal quantum number is in range
!
      IF (N<1 .OR. N>8) THEN
         WRITE (*, 300)
         WRITE (*, 301) N
         STOP
      ENDIF
!
!   Determine the azimuthal quantum number
!
      IF (KAPPA > 0) THEN
         L = KAPPA
      ELSE IF (KAPPA == 0) THEN
         WRITE (*, 300)
         WRITE (*, 302)
         STOP
      ELSE
         L = (-KAPPA) - 1
      ENDIF
!
!   Ensure that the azimuthal quantum number is in range
!
      IF (L > N - 1) THEN
         WRITE (*, 300)
         WRITE (*, 303) KAPPA, N
         STOP
      ENDIF
!
!   Find the appropriate entry in the table
!
      LOC = (N*N - N)/2 + L + 1
      BETHEL = BETHE(LOC)
!
!   Determine the quantity in square brackets in eq. (1) of
!   Klarsfeld and Maquet
!
      TERM = -BETHEL
!
      IF (KAPPA > 0) THEN
         TERM = TERM - C402/DBLE(L*(L + L + 1))
      ELSE
         TERM = TERM + C402/DBLE((L + 1)*(L + L + 1))
         IF (KAPPA == (-1)) THEN
            ZALFA = Z/C
            FACTOR = DLOG(ZALFA)
            FACTOR = -(FACTOR + FACTOR)
            TERM = TERM + FACTOR + C401
         ENDIF
      ENDIF
!
      FZALFA = OVLFAC*TERM
!
      RETURN
!
  300 FORMAT('KLAMAQ:')
  301 FORMAT(' Principal quantum number, ',1I2,&
         ', should be in the range  1--8.')
  302 FORMAT(' Kappa is  0 .')
  303 FORMAT(' Kappa, ',1I3,', is out of range for n, ',1I2,'.')
      RETURN
!
      END SUBROUTINE KLAMAQ
