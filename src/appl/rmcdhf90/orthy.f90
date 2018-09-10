!***********************************************************************
!
      SUBROUTINE ORTHY(NW, JP, LSORT) 
!
! nw     - input, total number of orbitals
! jp     - input, orbital jp, usually updated from  solve. It is
!          the one after IORDER
! lsort  - input, logical, .t. then sort according in decreasing
!          self-consistency order.
!
! The following quantities are not used outside this routine
!
! kfixed - output, the number of fixed orbitals of the same kappa
!          (specified by nak(jp)) as the orbital jp.
! ktotal - output, the number of all orbitals of the same kappa as the
!          orbital jp.
! kindx  - output, an index array containing positions of the orbitals
!          in the decreasing self-consistency order. It has the size of
!          ktotal.
!
! Orbitals are grouped as fixed, unfixed (non-correlation, correlation).
! .For i=1 to kfixed, kindx(i) has the same order as array iorder (though
!  this can be modified by requiring the smallest principal quantum number
!  come first).
! .For i=kfixed+1 to knon, kindx(i) is ordered as SCN/E**2.
! .For i=knon+1 to ktotal, kindx(i) is ordered as SCN.
! .When there is degeneracy (same SCN/E**2 or same SCN), iorder is used.
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:57:18   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNW
      USE DEF_C
      USE ORB_C, ONLY: nak, e
      USE ORBA_C, ONLY: iorder
      USE CORRE_C, ONLY: LCORRE 
      USE FIXD_C, ONLY: lfix
      USE INVT_C, ONLY: noinvt
      USE IOUNIT_C
      USE SCF_C, ONLY: scnsty
      USE WAVE_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rint_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NW 
      INTEGER , INTENT(IN) :: JP 
      LOGICAL , INTENT(IN) :: LSORT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1 
      INTEGER , DIMENSION(NNNW) :: KINDX 
      INTEGER :: KFIXED, NAKJ, KTOTAL, I, KI, J, KNON, ISHIFT, KJ, LRAW, L, &
         NAKL, MTP0, KRAW, K, MTP 
      REAL(DOUBLE) :: EPS, OVRLAP, DNORM, FACTOR 
      LOGICAL :: CHECK 
!-----------------------------------------------
!
! For inverse - taken from orthor.f which was not used anywhere
! Not used currently
!-----------------------------------------------------------------------
      CHECK = .NOT.NOINVT(JP) 
 
      EPS = 0.01D0*ACCY 
      KFIXED = 0 
      NAKJ = NAK(JP) 
 
!-----------------------------------------------------------------------
! Find orbitals of the same kappa as orbital jp.
!-----------------------------------------------------------------------
      KTOTAL = 0 
      DO I = 1, NW 
!         Orbital jp is supposed to be an ordered one, so should be k
         KI = IORDER(I) 
         IF (NAK(KI) /= NAKJ) CYCLE  
         KTOTAL = KTOTAL + 1 
         KINDX(KTOTAL) = KI 
      END DO 
 
!-----------------------------------------------------------------------
! Find fixed orbitals and separate them from the rest.
!-----------------------------------------------------------------------
 
!      Place the fixed orbitals at the begining of the array, no
!      sorting is done here. And find the number of non-correlation
!      orbitals.
 
      J = 0 
      KNON = 0 
      DO I = 1, KTOTAL 
         KI = KINDX(I) 
         IF (LFIX(KI)) THEN 
 
            J = J + 1 
 
             ! Shift array elements rather than exchanging
            KINDX(I:1+J:(-1)) = KINDX(I-1:J:(-1)) 
!            kindx(i) = kindx(j)
 
            KINDX(J) = KI 
 
         ELSE IF (.NOT.LCORRE(KI)) THEN 
            KNON = KNON + 1 
         ENDIF 
      END DO 
      KFIXED = J 
 
      IF (LSORT) THEN 
 
!-----------------------------------------------------------------------
! Separate correlation from non-correlation orbitals
!-----------------------------------------------------------------------
 
         J = KFIXED 
         DO I = KFIXED + 1, KTOTAL 
            KI = KINDX(I) 
            IF (LCORRE(KI)) CYCLE  
 
            J = J + 1 
 
                ! Shift array elements rather than exchanging
            KINDX(I:1+J:(-1)) = KINDX(I-1:J:(-1)) 
!               kindx(i) = kindx(j)
 
            KINDX(J) = KI 
         END DO 
 
!-----------------------------------------------------------------------
! Sort non-correlation and correlation orbitals separately
! This can be done indepently since knon is known
!-----------------------------------------------------------------------
 
!           non-correlation, using criteria scnsty/E^2
 
         DO I = KFIXED + 1, KFIXED + KNON 
            KI = KINDX(I) 
            DO J = I + 1, KFIXED + KNON 
               KJ = KINDX(J) 
               IF (SCNSTY(KJ)/E(KJ)**2 >= SCNSTY(KI)/E(KI)**2) CYCLE  
 
                   ! No need to do shifting, exchanging is fine
               KINDX(J) = KI 
               KI = KJ 
            END DO 
            KINDX(I) = KI 
         END DO 
 
!           correlation orbitals, using criteria scnsty
 
         DO I = KFIXED + KNON + 1, KTOTAL 
            KI = KINDX(I) 
            DO J = I + 1, KTOTAL 
               KJ = KINDX(J) 
               IF (SCNSTY(KJ) >= SCNSTY(KI)) CYCLE  
 
                   ! No need to do shifting, exchanging is fine
               KINDX(J) = KI 
               KI = KJ 
            END DO 
            KINDX(I) = KI 
         END DO 
 
      ENDIF 
 
! Finished sorting.
! Schmidt orthogonalize all orbitals of the same kappa
! The fixed orbitals are not changed
 
      DO LRAW = KFIXED + 1, KTOTAL 
         L = KINDX(LRAW) 
 
         NAKL = NAK(L) 
 
         MTP0 = MF(L) 
 
         DO KRAW = 1, LRAW - 1 
            K = KINDX(KRAW) 
            OVRLAP = RINT(L,K,0) 
 
!   Schmidt orthogonalise
 
            PZ(L) = PZ(L) - OVRLAP*PZ(K) 
            MTP = MAX(MF(L),MF(K)) 
            MTP0 = MAX(MTP0,MF(K)) 
 
            PF(:MTP,L) = PF(:MTP,L) - OVRLAP*PF(:MTP,K) 
            QF(:MTP,L) = QF(:MTP,L) - OVRLAP*QF(:MTP,K) 
         END DO 
 
!   Normalise
 
         MTP = MTP0 
 
         MF(L) = MTP 
         DNORM = RINT(L,L,0) 
         FACTOR = 1.D0/SQRT(DNORM) 
 
         IF (PZ(L) < 0.D0) FACTOR = -FACTOR 
 
         PZ(L) = FACTOR*PZ(L) 
         PF(2:MTP,L) = FACTOR*PF(2:MTP,L) 
         QF(2:MTP,L) = FACTOR*QF(2:MTP,L) 
 
!   Find new MF(L)
 
         MTP = MTP + 1 
   20    CONTINUE 
         MTP = MTP - 1 
         IF (ABS(PF(MTP,L)) < EPS) THEN 
            PF(MTP,L) = 0.D0 
            QF(MTP,L) = 0.D0 
            GO TO 20 
         ELSE 
            MF(L) = MTP 
         ENDIF 
 
      END DO 
 
      RETURN  
      END SUBROUTINE ORTHY 
