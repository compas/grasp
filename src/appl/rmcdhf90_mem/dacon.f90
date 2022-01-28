
!***********************************************************************
!                                                                      *
      SUBROUTINE DACON
!                                                                      *
!   This  routine  includes  the  contribution from the off-diagonal   *
!   I(a,b) integrals in the 'exchange' term.                           *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT.                                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE grid_C
      USE npot_C
      USE orb_C
      USE pote_C
      USE scf_C
      USE tatb_C
      USE wave_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IORB, MFI, I
      REAL(DOUBLE) :: TWOC, COEFF, FK, RPORII, PFI, QFI, ZBCI
!-----------------------------------------------
!
!
      TWOC = C + C
!
      DO K = 1, NDCOF
!
         IORB = NDA(K)
         CALL DPBDT (IORB)
         MFI = MF(IORB)
         COEFF = DA(K)
         FK = DBLE(NAK(IORB))
!
         DO I = 2, MFI
            RPORII = 1.0D0/(H*RPOR(I))
            PFI = PF(I,IORB)
            QFI = QF(I,IORB)
            ZBCI = ZZ(I)/C
            XP(I) = XP(I) + COEFF*(TA(I)*RPORII+FK*PFI-(TWOC*R(I)+ZBCI)*QFI)
            XQ(I) = XQ(I) + COEFF*(TB(I)*RPORII-FK*QFI+ZBCI*PFI)
         END DO
!
      END DO
!
      RETURN
      END SUBROUTINE DACON
