!***********************************************************************
!                                                                      *
      SUBROUTINE DAMPOR(J, INV, ODAMPJ)
!                                                                      *
!   This subroutine damps the orbital wave function with index J. it   *
!   also  stores  the  previous  determination  of this orbital.       *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:22:36   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP
      USE damp_C
      USE def_C
      USE grid_C
      USE int_C
      USE invt_C
      USE scf_C
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
      INTEGER, INTENT(INOUT) :: INV
      REAL(DOUBLE), INTENT(IN) :: ODAMPJ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTPO, MTPN, MTP, I, NNCFF, MFJ
      REAL(DOUBLE) :: EPS, FACTOR, POLDI, QOLDI, DNORM, DNFAC, SGN
      LOGICAL :: CHECK
!-----------------------------------------------
!
!
!   Initialization
!
!ww      EPS = 0.1D 00*ACCY
      EPS = 0.01D00*ACCY
      CHECK = .NOT.NOINVT(J)
!
!   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
!   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
!   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
!   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
!   dampck
!
!   Store previous determination of this orbital in arrays P, Q
!
!XHH odampj goes to the argument
!      ODAMPJ = ABS (ODAMP(J))
!
      IF (ODAMPJ > EPS) THEN
!
         FACTOR = 1.0D00 - ODAMPJ
!
         PZ(J) = FACTOR*P0 + ODAMPJ*PZ(J)
!
         MTPO = MF(J)
         MTPN = MTP0
         MTP0 = MTPO
!
         MTP = MAX(MTPN,MTPO)
         DO I = 1, MTP
            POLDI = PF(I,J)
            PF(I,J) = FACTOR*P(I) + ODAMPJ*PF(I,J)
            P(I) = POLDI
            QOLDI = QF(I,J)
            QF(I,J) = FACTOR*Q(I) + ODAMPJ*QF(I,J)
            Q(I) = QOLDI
         END DO
!
!   Compute normalization factor
!
         MF(J) = MTP
         DNORM = RINT(J,J,0)
         DNFAC = 1.0D00/SQRT(DNORM)
!
!   Determine if inversion is necessary
!
         IF (CHECK) THEN
            CALL COUNT (PF(:NNNP,J), MTP, NNCFF, SGN)
            IF (SGN < 0.0D00) THEN
               INV = INV + 1
               DNFAC = -DNFAC
            ENDIF
         ENDIF
!
         PZ(J) = PZ(J)*DNFAC
         PF(:MTP,J) = DNFAC*PF(:MTP,J)
         QF(:MTP,J) = DNFAC*QF(:MTP,J)
!
!   Find new MF(J)
!
         MFJ = MTP + 1
    3    CONTINUE
         MFJ = MFJ - 1
         IF (ABS(PF(MFJ,J)) < EPS) THEN
            PF(MFJ,J) = 0.0D00
            QF(MFJ,J) = 0.0D00
            GO TO 3
         ELSE
            MF(J) = MFJ
         ENDIF
!
      ELSE
!
         PZ(J) = P0
!
         MTPO = MF(J)
         MTPN = MTP0
         MTP0 = MTPO
!
         MTP = MAX(MTPN,MTPO)
         DO I = 1, MTP
            POLDI = PF(I,J)
            PF(I,J) = P(I)
            P(I) = POLDI
            QOLDI = QF(I,J)
            QF(I,J) = Q(I)
            Q(I) = QOLDI
         END DO
!
         MF(J) = MTP
!
      ENDIF
!
      RETURN
      END SUBROUTINE DAMPOR
