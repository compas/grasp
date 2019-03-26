!***********************************************************************
!                                                                      *
      SUBROUTINE HMOUT(IMCDF)
!                                                                      *
!   Routine for printing the Hamiltonian matrix.  File IMCDF must be   *
!   positioned correctly before a call is made to this module.         *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE hmat_C
      USE orb_C,           ONLY: ncf, nw, iqa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: IMCDF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IC, NELC, IR, LENTHR, LENTHC
      REAL(DOUBLE) :: ELSTO
      CHARACTER    :: CIR*8, CIC*8
!-----------------------------------------------
!
      DO IC = 1, NCF
         READ (IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), (IROW(IR),IR=1,NELC)
         DO IR = 1, NELC
            CALL CONVRT (IROW(IR), CIR, LENTHR)
            CALL CONVRT (IC, CIC, LENTHC)
            WRITE (99, 300) CIR(1:LENTHR), CIC(1:LENTHC), EMT(IR)
         END DO
      END DO
!
  300 FORMAT(' H (',A,',',A,') = ',1P,1D19.12)
      RETURN
!
      END SUBROUTINE HMOUT
