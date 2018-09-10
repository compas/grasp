!
!     ------------------------------------------------------------------
!     P A M T M T
!     ------------------------------------------------------------------
!
      SUBROUTINE PAMTMT(X, T, WORK, NORB) 
!
! GENERATE PER AKE'S T MATRIX FROM A
! ORBITAL ROTATION MATRIX X
!
! T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
! MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
!
!         TL = 1 - L
!         TU = U ** -1
!
! WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
! X :
!         X = L * U
! WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
! DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
!
! JEPPE OLSEN OCTOBER 1988
!
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lulu_I 
      USE setvec_I 
      USE wrtmat_I 
      USE invmat_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NORB 
      REAL(DOUBLE)  :: X(NORB,NORB) 
      REAL(DOUBLE)  :: T(NORB,NORB) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTEST, KLFREE, KLL, KLU, I, J 
!-----------------------------------------------
! DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
!
      NTEST = 0 
!. Allocate local memory
      KLFREE = 1 
!     KLL = KFLREE
      KLL = KLFREE 
      KLFREE = KLL + NORB*(NORB + 1)/2 
      KLU = KLFREE 
      KLFREE = KLU + NORB**2 
!.LU factorize X
      CALL LULU (X, WORK(KLL), WORK(KLU), NORB) 
!.Expand U to full matrix
      CALL SETVEC (T, 0.0D0, NORB**2) 
      DO I = 1, NORB 
         DO J = I, NORB 
            T(I,J) = WORK(KLU-1+J*(J-1)/2+I) 
         END DO 
      END DO 
      IF (NTEST >= 10) THEN 
         WRITE (6, *) ' MATRIX TO BE INVERTED ' 
         CALL WRTMAT (T, NORB, NORB, NORB, NORB) 
      ENDIF 
!.Invert U
      CALL INVMAT (T, WORK(KLU), NORB, NORB) 
      IF (NTEST >= 10) THEN 
         WRITE (6, *) ' INVERTED MATRIX ' 
         CALL WRTMAT (T, NORB, NORB, NORB, NORB) 
      ENDIF 
!.Subtract L
      DO I = 1, NORB 
         T(I,:I-1) = -WORK(KLL+I*(I-1)/2:I-2+KLL+I*(I-1)/2) 
      END DO 
!
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' INPUT X MATRIX ' 
         CALL WRTMAT (X, NORB, NORB, NORB, NORB) 
         WRITE (6, *) ' T MATRIX ' 
         CALL WRTMAT (T, NORB, NORB, NORB, NORB) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE PAMTMT 
