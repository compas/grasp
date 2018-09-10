!
!     ------------------------------------------------------------------
!     U L L A
!     ------------------------------------------------------------------
!
      SUBROUTINE ULLA(A, U, L, NDIM, SCR) 
!
! Obtain U L decomposition of matrix A
!   A = U L
!
! Note that L and U are returned in full matrix form
!
! Quick and dirty routine, Jeppe Olsen, November 1991
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
      USE lulu_I 
      USE setvec_I 
      USE wrtmat_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NDIM 
      REAL(DOUBLE), INTENT(IN) :: A(NDIM,NDIM) 
      REAL(DOUBLE)  :: U(NDIM,NDIM) 
      REAL(DOUBLE)  :: L(NDIM,NDIM) 
      REAL(DOUBLE)  :: SCR(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KLFREE, KLPAP, I, J, KLL, KLU, IMAX, IMIN, NTEST 
!-----------------------------------------------
!
!
! In order to change into LU form introduce the orthogonal matrix
!    p(i,j) = delta(i,ndim-i)
! and rewrite
! A = P P A P P  = P L U P = PLP PUP
! where LU is an LU decomposition of PAP, since PLP is upper
! tringular AND PUP is lower traingular we have obtained the goal
!
! 1 : PAP in scr(klPAP)
!
      KLFREE = 1 
      KLPAP = KLFREE 
      KLFREE = KLFREE + NDIM**2 
!
      DO I = 1, NDIM 
         SCR(KLPAP-1+I:NDIM*(NDIM-1)+KLPAP-1+I:NDIM) = A(NDIM+1-I,NDIM:1:(-1)) 
      END DO 
! 2 : Lu decompose PAP
      KLL = KLFREE 
      KLFREE = KLFREE + NDIM*(NDIM + 1)/2 
!
      KLU = KLFREE 
      KLFREE = KLFREE + NDIM*(NDIM + 1)/2 
      CALL LULU (SCR(KLPAP), SCR(KLL), SCR(KLU), NDIM) 
!          LULU(A,L,U,NDIM)
! Storage modes
!   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
!   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
!
!. 3 : Obtain U as PLP and L as PUP
      CALL SETVEC (U, 0.0D0, NDIM**2) 
      CALL SETVEC (L, 0.0D0, NDIM**2) 
!
      DO IMAX = 1, NDIM 
         DO IMIN = 1, IMAX 
            U(IMIN,IMAX) = SCR(KLL-1+(NDIM+1-IMIN)*(NDIM+1-IMIN-1)/2+(NDIM+1-&
               IMAX)) 
            L(IMAX,IMIN) = SCR(KLU-1+(NDIM+1-IMIN)*(NDIM+1-IMIN-1)/2+(NDIM+1-&
               IMAX)) 
         END DO 
      END DO 
!
      NTEST = 0 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' U and L from Ulla ' 
         CALL WRTMAT (U, NDIM, NDIM, NDIM, NDIM) 
         CALL WRTMAT (L, NDIM, NDIM, NDIM, NDIM) 
      ENDIF 
!
 
      RETURN  
      END SUBROUTINE ULLA 
