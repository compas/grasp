!***********************************************************************
!                                                                      *
      SUBROUTINE genintbreit2 (myid, nprocs, NB, j2max)
!                                                                      *
!  Input:
!     myid, nprocs
!  Output:
!   N - Number of integrals
!   j2max - max of 2*J
!
!       Generate the list of Breit Integrals of type 2                 *
!       that could arise from a set of orbitals                        *
!       of orbitals.                                                   *
!                                                                      *
!       This routine is similar to genintrk                            *
!                                                                      *
!     Written by Per Jonsson            October 2014                   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer 
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE bilst_C
      USE orb_C, ONLY: nkj, nw
      USE kkstartbreit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE brintf_I
      USE triangbreit2_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: myid, nprocs
      INTEGER             :: NB, j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: GEN
      INTEGER :: I, KEY, K, IA, IB, IC, ID
!-----------------------------------------------
!

      KEY = NW + 1
      KSTARTBREIT2(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintbreit2: KMAX too small'
      ENDIF
!
!   Make the breit integrals. Note that there are no symmetry relations that
!   can be utilized. See paper by Grasnt.
!   When GEN is false, sweep through to find dimension
!
      GEN = .FALSE.

  999 NB = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
            DO IB = 1, NW
               DO IC = 1, NW
                  DO ID = 1, NW
                     IF (TRIANGBREIT2(IA,IB,IC,ID,K)) THEN
                        NB = NB + 1
                        IF (GEN .AND.(MOD(NB,nprocs) .EQ. myid)) THEN
                           INDTP2(NB) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                           VALTP2(NB) = BRINTF (2,IA,IB,IC,ID,K)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         KSTARTBREIT2(K+1) = NB + 1
      ENDDO
!
!     Allocate memory for integral book keeping
!
      IF (.NOT. GEN) THEN
         CALL ALLOC (INDTP2,NB,'INDT2','genintbreit2')
         CALL ALLOC (VALTP2,NB,'VALT2','genintbreit2')

! Initialization is necessary in the mpi version

         DO i = 1, NB
            INDTP2(i) = 0
            VALTP2(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0)                                        &
            PRINT *, 'Computing',NB,' Breit integrals of type 2'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END SUBROUTINE genintbreit2
