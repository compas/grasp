!***********************************************************************
!                                                                      *
      SUBROUTINE ORTHOR(J, INV) 
!
!                                                                      *
!   This routine Schmidt orthogonalizes orbital  J  to  all orbitals   *
!   which  have  better  self-consistency.  Note that fixed orbitals   *
!   have the best self-consistency.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
! INV not initialized here but used as
! INV = INV + 1
! A bug regarding IORDER fixed and new schemes introduced. See comments
! Anyway this routine was not used anywhere.
!XHH 1997.02.21
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:54:29   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNP
      USE DEF_C 
      USE GRID_C
      USE OVL_C 
      USE ORB_C 
      USE ORBA_C
      USE INVT_C
      USE SCF_C
      USE IOUNIT_C
      USE FIXD_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rint_I 
!      USE count_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J 
      INTEGER , INTENT(INOUT) :: INV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KRAW, K, MTP, I, NNCFF, JGG
      REAL(DOUBLE) :: EPS, OVRLAP, DNORM, FACTOR, SGN 
      LOGICAL :: CHECK, CHANGED 
!-----------------------------------------------
!
!XHH Added common /orba/ and /iounit/
 
!
!ww      EPS = ACCY*0.1D 00
      EPS = ACCY*0.01D00 
!
      CHECK = .NOT.NOINVT(J) 
!
!XHH
! Bug fixed. Here J is actually IORDER(J_raw), thus K should be
! treated the same way.
      CHANGED = .FALSE. 
!      DO 4 K = 1,NW
      DO KRAW = 1, NW 
         K = IORDER(KRAW) 
!         write(istde,*) '***',kraw,k,np(k),nh(k),scnsty(k),'***'
!
!XHH orbitals with higher self-consistency are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) .AND.
!     :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
!XHH All orbitals are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) ) THEN
!XHH orbitals before the current and unchanged ones are considered
         IF (.NOT.(NAK(K)==NAK(J) .AND. (K<J .OR. LFIX(K)))) CYCLE  
 
         CHANGED = .TRUE. 
!
!   Compute overlap
!
         OVRLAP = RINT(J,K,0) 
!
!   Schmidt orthogonalise
!
         PZ(J) = PZ(J) - OVRLAP*PZ(K)
         MTP = MAX(MF(J),MF(K)) 
         PF(:MTP,J) = PF(:MTP,J) - OVRLAP*PF(:MTP,K) 
         QF(:MTP,J) = QF(:MTP,J) - OVRLAP*QF(:MTP,K) 
      END DO 
!XHH
      IF (CHANGED) THEN 
!
!   Normalise
!
         MF(J) = MTP 
         DNORM = RINT(J,J,0) 
         FACTOR = 1.0D00/SQRT(DNORM) 
!
!   Determine if inversion is necessary
!
         IF (CHECK) THEN 
            CALL COUNT (PF(:NNNP,J), MTP, NNCFF, SGN) 
            IF (SGN < 0.0D00) THEN 
               INV = INV + 1 
               FACTOR = -FACTOR 
            ENDIF 
         ENDIF 
!
!   Perform normalization and/or inversion
!
         PZ(J) = FACTOR*PZ(J)
         PF(2:MTP,J) = FACTOR*PF(2:MTP,J) 
         QF(2:MTP,J) = FACTOR*QF(2:MTP,J) 
!
!   Find new MF(J)
!
         MTP = MTP + 1 
    3    CONTINUE 
         MTP = MTP - 1 
         IF (ABS(PF(MTP,J)) < EPS) THEN 
            PF(MTP,J) = 0.0D00 
            QF(MTP,J) = 0.0D00 
            GO TO 3 
         ELSE 
            MF(J) = MTP 
         ENDIF 
!
!XHH two lines moved forward and added an ENDIF
!         ENDIF
!    4 CONTINUE
      ENDIF 
!
      RETURN  
      END SUBROUTINE ORTHOR 
