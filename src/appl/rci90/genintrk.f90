!***********************************************************************
!                                                                      *
      SUBROUTINE genintrk (myid, nprocs, N, j2max)
!                                                                      *
!  Input:
!     myid, nprocs
!  Output:
!   N - Number of integrals
!   j2max - max of 2*J
!
!       Generate the list of Rk Integrals that could arise from a set  *
!       of orbitals.                                                   *
!                                                                      *
!     Written by Per Jonsson                                           *
!   MPI version by Xinghong He            Last revision: 22 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE cteilsrk_C
      USE kkstart_C
      USE orb_C
      USE jlabl_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE slater_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: myid, nprocs
      INTEGER, INTENT(OUT) :: N, j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: GEN,TRIANGRK
      INTEGER :: key, i, k, ia, ib, ic, id
!-----------------------------------------------------------------------

      KEY = NW + 1
      KSTART(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintrk: KMAX too small'
      ENDIF
!
!   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
!   When GEN is false, sweep through to find dimension
!
      GEN = .FALSE.

  999 N = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
            DO IB = IA, NW
               DO IC = IA, NW
                  IF (TRIANGRK(NKL(IA),K,NKL(IC))) THEN
                     DO ID = IB, NW
                        IF (TRIANGRK(NKL(IB),K,NKL(ID))) THEN
                           N = N + 1
                           IF (GEN .AND. (MOD(N,nprocs) .EQ. myid)) THEN
                              INDTEIRK(N) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                              VALTEIRK(N) = SLATER (IA,IB,IC,ID,K)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         KSTART(K+1) = N + 1
      ENDDO
!
!     Allocate memory for integral book keeping
!
      IF (.NOT. GEN) THEN
         CALL ALLOC (INDTEIRK,N,'INDTEIRK', 'GENINTRK')
         CALL ALLOC (VALTEIRK,N,'VALTEIRK', 'GENINTRK')

! Initialization is necessary in the mpi version

         DO i = 1, N
            INDTEIRK(i) = 0
            VALTEIRK(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0)                                               &
     &      PRINT *, 'Allocating space for ',N,' Rk integrals'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END
