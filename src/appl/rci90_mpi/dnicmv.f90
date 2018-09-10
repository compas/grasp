!***********************************************************************
!                                                                      *
      SUBROUTINE DNICMV(N, M, B, C) 
!                                                                      *
!   Matrix-matrix product: C = AB.  The lower triangle of the  (NxN)   *
!   matrix is assumed available in packed form in the array EMT. The   *
!   matrices B and C are (NxM).                                        *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  TRSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
!               [BLAS]: DAXPY/SAXPY, DDOT/SDOT.                        *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 09 Oct 1992   *
!   Block version by Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE hmat_C,          ONLY: EMT
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dinit_I 
!      USE dmerge_dnicmv_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N 
      INTEGER, INTENT(IN) :: M 
      REAL(DOUBLE), DIMENSION(N,M) :: B
      REAL(DOUBLE), DIMENSION(N,M) :: C
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IBEG, IEND, ICOL, NELC, IV 
      REAL(DOUBLE) :: DIAG, DL 
!-----------------------------------------------
 
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- there is no overdimensioning
 
      CALL DINIT (N*M, 0.0D00, C, 1) 
 
      IBEG = 1 
      IEND = 0 
      DO ICOL = MYID + 1, N, NPROCS 
         IEND = IEND + ICOL 
         NELC = IEND - IBEG + 1 
         DO IV = 1, M 
            DIAG = C(ICOL,IV) + EMT(IEND)*B(ICOL,IV) 
            CALL DMERGE_DNICMV (NELC - 1, B(1:N,IV), C(1:N,IV), &
                                EMT(IBEG:IEND), B(ICOL,IV), DL) 
            C(ICOL,IV) = DIAG + DL 
         END DO 
         IBEG = IEND + 1 
      END DO 

      CALL gdsummpi (C,N*M)
 
      RETURN  
      END SUBROUTINE DNICMV 


 
!***********************************************************************
!                                                                      *
      SUBROUTINE DMERGE_DNICMV(N, DB, DC, DA, DCONST, DL) 
!
!  Used by dnimcv
!  Developed from dmerge. The only diff is: idy not needed here
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N 
      REAL(DOUBLE), INTENT(IN) :: DCONST 
      REAL(DOUBLE), INTENT(OUT) :: DL 
      REAL(DOUBLE), DIMENSION(*), INTENT(IN) :: DB
      REAL(DOUBLE), DIMENSION(*), INTENT(INOUT) :: DC
      REAL(DOUBLE), DIMENSION(N), INTENT(IN) :: DA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: DSUM 
!-----------------------------------------------
 
      DSUM = 0.D0 
      DSUM = DOT_PRODUCT(DA,DB(:N)) 
      DC(:N) = DC(:N) + DCONST*DA 
      DL = DSUM 
 
      RETURN  
      END SUBROUTINE DMERGE_DNICMV 
