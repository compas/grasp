!***********************************************************************
      SUBROUTINE SPICMVmpi (N,M,B,C)
!     IMPLICIT REAL*8          (A-H,O-Z)

!   This routine now works for both rscfmpivu and rcimpivu. By removing
!   the include statement and the call to gdsummpi and setting myid=0,
!   nprocs=1, it also works for the corresponding serial program.
!   Matrix is stored in the mode of upper-triangle-by columns, or
!   you can say lower-triangle-by-rows. 98-08-06
!
!   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
!   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
!   COMMON/HMAT/.                                                      *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT                                 *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 19 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jul 1998   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE mpi_C
      USE HMAT_C,          ONLY: EMT, IENDC, IROW, NELMNT
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dinit_I
!      USE dmerge_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: M
      REAL(DOUBLE) :: B(N,M)
      REAL(DOUBLE) :: C(N,M)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICOL, IBEG, IEND, NELC, IV
      REAL(DOUBLE) :: DIAG, DL
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- no overdimensioning
!
      CALL DINIT (N*M, 0.0D00, C, 1)

      ibeg = 1
      DO ICOL = myid + 1, N, nprocs
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         IEND = IENDC(ICOL)
         NELC = IEND - IBEG + 1
         DO IV = 1, M
            DIAG =  C(ICOL,IV) + EMT(iend)*B(icol,IV)
            CALL DMERGE (NELC-1,B(:N,IV),C(:N,IV),               &
                         IROW(IBEG),EMT(IBEG),B(ICOL,IV),DL)
            C(ICOL,IV) = DIAG + DL
         ENDDO
         ibeg = iend + 1
      ENDDO

      CALL gdsummpi (C, N*M)

      RETURN
      END SUBROUTINE SPICMVmpi
