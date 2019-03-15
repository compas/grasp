!***********************************************************************
!
      SUBROUTINE GETALD
!
!   Interactively determines the data governing AL problem.
!
!   Quantities to be determined:
!     ncmin, wt(1:ncf), ucf(1:nw), nscf, nsic, orthst
!
!   Call(s) to: [LIB92]: ALLOC, IQ.
!
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992
!   Block version by Xinghong He          Last revision: 13 Jul 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE def_C
      USE fixd_C
      USE orb_C
      USE orthct_C
      USE scf_C
      USE iounit_C
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getaldwt_I
      USE iq_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM, J, I
      REAL(DOUBLE) :: SUM
      LOGICAL :: GETYN, YES
!-----------------------------------------------
!
      WRITE (ISTDE, *) '(E)AL type calculation; H(DC) will not be ', &
         'diagonalised;'
      WRITE (ISTDE, *) 'getald ...'
      WRITE (ISTDE, *) 'ncf=', NCF

      CALL ALLOC (WT, NCF, 'WT', 'GETALD')

      CALL GETALDWT (NCF, WT)

      DO J = 1, NW
         SUM = 0.D0
         DO I = 1, NCF
            SUM = SUM + WT(I)*DBLE(IQ(J,I))
         END DO
         UCF(J) = SUM
      END DO

      NCMIN = 0
      NSCF = 12
      NSIC = 2 + (NW - NFIX)/4
      ORTHST = .FALSE.

      RETURN
      END SUBROUTINE GETALD
