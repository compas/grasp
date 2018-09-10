!
!
!     ------------------------------------------------------------------
!     L U L U
!     ------------------------------------------------------------------
!
      SUBROUTINE LULU(A, L, U, NDIM) 
!
! LU DECOMPOSITION OF MATRIX A
!
!     A = L * U
!
! WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
! UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
!
! L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
!
!   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
!
!   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
!
! THIS ADRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
! FOR L AND  OVER ROWS FOR U .
!
!
! NO PIVOTING IS DONE HERE , SO THE SCHEME GOES :
!
!     LOOP OVER R=1, NDIM
!        LOOP OVER J = R, NDIM
!          U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
!        END OF LOOP OVER J
!
!        LOOP OVER I = R+1, NDIM
!          L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
!        END OF LOOP OVER I
!     END OF LOOP OVER R
!
! JEPPE OLSEN , OCTOBER 1988
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE inprod_I 
      USE prsym_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NDIM 
      REAL(DOUBLE), INTENT(IN) :: A(NDIM,NDIM) 
      REAL(DOUBLE)  :: L(*) 
      REAL(DOUBLE)  :: U(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: R, J, I, NTEST 
      REAL(DOUBLE) :: XFACI 
!-----------------------------------------------
!
      DO R = 1, NDIM 
!
         DO J = R, NDIM 
            U(J*(J-1)/2+R) = A(R,J) - INPROD(L(R*(R-1)/2+1),U(J*(J-1)/2+1),R-1) 
         END DO 
!
         XFACI = 1.0D0/U(R*(R+1)/2) 
         L(R*(R+1)/2) = 1.0D0 
         DO I = R + 1, NDIM 
            L(I*(I-1)/2+R) = (A(I,R)-INPROD(L(I*(I-1)/2+1),U(R*(R-1)/2+1),R-1))&
               *XFACI 
         END DO 
!
      END DO 
!
      NTEST = 0 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' L MATRIX ' 
         CALL PRSYM (L, NDIM) 
         WRITE (6, *) ' U MATRIX ( TRANSPOSED ) ' 
         CALL PRSYM (U, NDIM) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE LULU 
