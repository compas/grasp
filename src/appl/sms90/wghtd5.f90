!***********************************************************************
!                                                                      *
      SUBROUTINE WGHTD5
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: ALLOC, DALLOC, ISPAR, ITJPO.                           *
!                                                                      *
!                                          Last updated: 02 Nov 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  08:39:50   2/21/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE EIGV_C
      USE JLABL_C
      USE ORB_C
      USE PRNT_C
      USE SYMA_C
      USE itjpo_I
      use ispar_I
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(5) :: ICONF
      REAL(DOUBLE),dimension(5)::wght
      INTEGER :: NELT, NVEX, IV, ICF, IFIRST, I, M, L, IP
      integer, dimension(:), pointer :: NEXT;
      REAL(DOUBLE),dimension(:), pointer :: WT
!-----------------------------------------------
!
!   Allocate storage for local arrays
!
      CALL ALLOC (WT, NCF, 'WT', 'WGHTD5' )
      CALL ALLOC (NEXT, NCF, 'NEXT', 'WGHTD5')
!
      WRITE (24, 300)
!
      IF (NCF < 5) NELT = NCF
!
      NVEX = NVEC
!
      DO IV = 1, NVEX
!
         ICF = IVEC(IV)
!
!   Set up linked list of weights
!
         NEXT(1) = 0
         WT(1) = EVEC(1 + (IV - 1)*NCF)**2
         IFIRST = 1
         L4: DO I = 2, NCF
            M = IFIRST
            L = 0
            WT(I) = EVEC(I + (IV - 1)*NCF)**2
            IF (WT(I) > WT(M)) THEN
               IF (L /= 0) GO TO 3
               NEXT(I) = IFIRST
               IFIRST = I
               CYCLE  L4
            ENDIF
            L = M
            M = NEXT(L)
            DO WHILE(M /= 0)
               IF (WT(I) > WT(M)) THEN
                  IF (L /= 0) GO TO 3
                  NEXT(I) = IFIRST
                  IFIRST = I
                  CYCLE  L4
               ENDIF
               L = M
               M = NEXT(L)
            END DO
    3       CONTINUE
            NEXT(I) = NEXT(L)
            NEXT(L) = I
         END DO L4
!
!   Print first five elements of list.
!
         M = IFIRST
         I = 0
         IF (ITJPO(M)==IATJPO(IV) .AND. ISPAR(M)==IASPAR(IV)) THEN
            I = I + 1
            WGHT(I) = WT(M)
            ICONF(I) = M
         ENDIF
         M = NEXT(M)
         DO WHILE(M/=0 .AND. I<5)
            IF (ITJPO(M)==IATJPO(IV) .AND. ISPAR(M)==IASPAR(IV)) THEN
               I = I + 1
               WGHT(I) = WT(M)
               ICONF(I) = M
            ENDIF
            M = NEXT(M)
         END DO
         IP = (IASPAR(IV) + 3)/2
         NELT = MIN(I,5)
         WRITE (24, 301) ICF, JLBL(IATJPO(IV)), JLBP(IP), (WGHT(I),I=1,NELT)
         WRITE (24, 302) (ICONF(I),I=1,NELT)
      END DO
!
!   Deallocate storage for local arrays
!
      CALL DALLOC (WT, 'WT', 'WGHTD5')
      CALL DALLOC (NEXT, 'NEXT', 'WGHTD5')
!
      RETURN
!
  300 FORMAT(/,'Weights of major contributors to ASF:'/,/,&
         'Level J Parity      CSF contributions'/)
  301 FORMAT(I3,2X,2A4,5(D12.4))
  302 FORMAT(13X,5I12)
      RETURN
!
      END SUBROUTINE WGHTD5
