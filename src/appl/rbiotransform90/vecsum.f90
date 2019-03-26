!
!     ------------------------------------------------------------------
!     V E C S U M
!     ------------------------------------------------------------------
!
      SUBROUTINE VECSUM(C, A, B, FACA, FACB, NDIM)
!
!     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NDIM
      REAL(DOUBLE), INTENT(IN) :: FACA
      REAL(DOUBLE), INTENT(IN) :: FACB
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(IN)  :: A
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(IN)  :: B
      REAL(DOUBLE), DIMENSION(NDIM), INTENT(OUT) :: C
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: S
!-----------------------------------------------
!
      IF (FACA/=0.0D0 .AND. FACB/=0.0D0) THEN
         C(:NDIM) = FACA*A(:NDIM) + FACB*B(:NDIM)
!
      ELSE IF (FACA==0.0D0 .AND. FACB/=0.0D0) THEN
         C(:NDIM) = FACB*B(:NDIM)
!
      ELSE IF (FACA/=0.0D0 .AND. FACB==0.0D0) THEN
         C(:NDIM) = FACA*A(:NDIM)
!
      ELSE IF (FACA==0.0D0 .AND. FACB==0.0D0) THEN
         C(:NDIM) = 0.0D0
      ENDIF
!
      RETURN
      END SUBROUTINE VECSUM
