!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION FZALF (N, KAPPA, Z)
!                                                                      *
!   An estimate of the function  F (Z*\alpha) is computed here.        *
!                                                                      *
!   Call(s) to: [RCI92]: KLAMAQ, MOHR.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1990   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mohr_I
      USE klamaq_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: KAPPA
      REAL(DOUBLE)  :: Z
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEFF
      REAL(DOUBLE) :: VALUE
!-----------------------------------------------
!
      IF (N <= 2) THEN
         CALL MOHR (N, KAPPA, Z, VALUE)
      ELSE
         IF (KAPPA==(-1) .OR. KAPPA==1 .OR. KAPPA==(-2)) THEN
            NEFF = 2
            CALL MOHR (NEFF, KAPPA, Z, VALUE)
         ELSE
            CALL KLAMAQ (N, KAPPA, Z, VALUE)
         ENDIF
      ENDIF
!
      FZALF = VALUE
!
      RETURN
      END FUNCTION FZALF
