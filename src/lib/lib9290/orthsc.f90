!***********************************************************************
!                                                                      *
      SUBROUTINE ORTHSC
!                                                                      *
!   This routine Schmidt orthogonalises radial wavefunctions.          *
!                                                                      *
!   Call(s) to: [LIB92]: RINT.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
!                                                                      *
! Normalization of the orbitals moved out of the inner loop
! XHH 1997.02.14
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:41:29   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW
      USE DEBUG_C
      USE DEF_C
      USE GRID_C
      USE ORB_C
      USE WAVE_C, ONLY: PZ, MF, PF, QF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rint_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NNNW) :: J
      INTEGER :: L, NAKL, KOUNT, MTP0, K, MTP, I
      REAL(DOUBLE), DIMENSION(NNNW) :: OVLAP
      REAL(DOUBLE) :: EPS, OVRLAP, DNORM, FACTOR
      LOGICAL :: CHANGED
!-----------------------------------------------
!
!XHH
!
!   Set tabulated values of the radial wavefunction to zero
!   if they are less than EPS
!
      EPS = 0.01D00*ACCY
!
!   Determine the number of interesting overlaps
!
!XHH wasting time, no use at all
!      NOVL = 0
!      DO 2 K = 1,NW-1
!         NAKK = NAK(K)
!         DO 1 L = K+1,NW
!            IF (NAKK .EQ. NAK(L)) NOVL = NOVL+1
!    1    CONTINUE
!    2 CONTINUE
!
!      IF (NOVL .EQ. 0) RETURN
!
      DO L = 2, NW
!
         NAKL = NAK(L)
         KOUNT = 0
!XHH MTP0 introduced to count the maximum number of points during
!    orthogonalization of the L-th orbital to other orbitals
!    The logical variable changed is initialized to .F.

         MTP0 = MF(L)
         CHANGED = .FALSE.
!
         DO K = 1, L - 1
!
            IF (NAK(K) /= NAKL) CYCLE
!XHH
            CHANGED = .TRUE.
!
!   Compute overlap
!
            OVRLAP = RINT(L,K,0)
!
!   Schmidt orthogonalise
!
            KOUNT = KOUNT + 1
            J(KOUNT) = K
            OVLAP(KOUNT) = OVRLAP
!
            PZ(L) = PZ(L) - OVRLAP*PZ(K)
            MTP = MAX(MF(L),MF(K))
            MTP0 = MAX(MTP0,MF(K))

            PF(:MTP,L) = PF(:MTP,L) - OVRLAP*PF(:MTP,K)
            QF(:MTP,L) = QF(:MTP,L) - OVRLAP*QF(:MTP,K)
         END DO
!
!   Normalise
!
!XHH Use MTP0 to replace MTP and only when the orbital is changed.
!    This is in accordance with the original version which had the
!    normalization etc within the inner K loop.

         IF (CHANGED) THEN
            MTP = MTP0

            MF(L) = MTP
            DNORM = RINT(L,L,0)
            FACTOR = 1.0D00/SQRT(DNORM)
!
            PZ(L) = FACTOR*PZ(L)
            PF(2:MTP,L) = FACTOR*PF(2:MTP,L)
            QF(2:MTP,L) = FACTOR*QF(2:MTP,L)
!
!   Find new MF(L)
!
            MTP = MTP + 1
    7       CONTINUE
            MTP = MTP - 1
!cjb       print *, '    MTP    = ',    MTP
!          print *, '        L  = ',        L
!          print *, ' PF(MTP,L) = ', PF(MTP,L)
!cjb Subscript out of range for array pf
           if ( MTP .GE. 1 ) then
!cjb
            IF (ABS(PF(MTP,L)) < EPS) THEN
               PF(MTP,L) = 0.0D00
               QF(MTP,L) = 0.0D00
               GO TO 7
            ELSE
               MF(L) = MTP
            ENDIF
!cjb endif
           endif
         ENDIF
         IF (.NOT.(LDBPR(3) .AND. KOUNT>0)) CYCLE

!XHH Moved ahead
!            ENDIF
!    8    CONTINUE
!
!   Print overlap information
!
!---------------------------------------------------------------
! debug 97.02.14 done
! Check the orthonormalisation
! KOUNT and OVLAP are destroyed
! should be commented out after seeing the values
!
!      KOUNT = 0
!      DO k = 1, L
!            IF (NAK(K) .EQ. NAKL) THEN
!               KOUNT = KOUNT+1
!               OVLAP(KOUNT) = RINT (L,K,0)
!               J(KOUNT) = K
!            ENDIF
!      ENDDO
!      WRITE (*,301)
!     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)
!---------------------------------------------------------------
         WRITE (99, 301) (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I=1,KOUNT)
!
      END DO
!
      RETURN
!
  301 FORMAT(1P,5(2X,1D10.3,' = <',1I2,1A2,'|',1I2,1A2,'>'))
      RETURN
!
      END SUBROUTINE ORTHSC
