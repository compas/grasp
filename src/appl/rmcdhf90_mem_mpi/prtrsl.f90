!***********************************************************************
!                                                                      *
      SUBROUTINE PRTRSL
!                                                                      *
!   This subroutine is now called only in GETOLD to print out all
!   orbitals so that the user can do cut-and-paste in supplying the
!   orbitals to be varied. I.e., LFIX(LOC) will always be FALSE and
!   are thus commented out.
! Xinghong HE 6 Apr 1998
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE orb_C
      USE fixd_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IEND, I, LOC, IBEG, LENTH
      CHARACTER :: RECORD*80, CNUM*2
!-----------------------------------------------
!
!
!   Print the list of all subshell radial wavefunctions
!
      IEND = 0
      DO I = 1, NW
!         LOC = IORDER(I)
         LOC = I
! For the commenting-out of the IF...ENDIF, see comments in the header
!         IF (.NOT. LFIX(LOC)) THEN
         IF (IEND > 75) THEN
            WRITE (ISTDE, *) RECORD(1:IEND)
            IEND = 0
         ENDIF
         IF (IEND > 0) THEN
            IBEG = IEND + 1
            IEND = IBEG
            RECORD(IBEG:IEND) = ' '
         ENDIF
         IBEG = IEND + 1
         CALL CONVRT (NP(LOC), CNUM, LENTH)
         IEND = IBEG + LENTH - 1
         RECORD(IBEG:IEND) = CNUM(1:LENTH)
         IF (NAK(LOC) < 0) THEN
            LENTH = 1
         ELSE
            LENTH = 2
         ENDIF
         IBEG = IEND + 1
         IEND = IBEG + LENTH - 1
         RECORD(IBEG:IEND) = NH(LOC)(1:LENTH)
!         ENDIF
      END DO
      IF (IEND > 1) WRITE (ISTDE, *) RECORD(1:IEND)
!
      RETURN
      END SUBROUTINE PRTRSL
