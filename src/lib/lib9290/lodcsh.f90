!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSH(NFILE, NCORE)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]:  PRSRCN, PRSRSL       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
!   Modified by C. F. Fischer to read only the header information
!
! Input:
!   nfile
!
! Output:
!   ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:32:39   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE DEF_C
      USE ORB_C
      USE TERMS_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE prsrsl_I
      USE prsrcn_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NFILE
      INTEGER  :: NCORE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NNNW) :: IOCC
      INTEGER :: IQADUM, NPEEL, I, J, NPJ, NAKJ, IOS, IERR
      CHARACTER :: STR*256
!-----------------------------------------------
!
!   Entry message
!
      WRITE (6, *) 'Loading CSF file ... Header only'
!
!   Get the list of subshells
!
      NW = 0
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      CALL PRSRSL (NFILE, 1)
      NCORE = NW
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      READ (NFILE, *)
      CALL PRSRSL (NFILE, 2)
      NPEEL = NW - NCORE
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      DO J = NCORE + 1, NW
         NPJ = NP(J)
         NAKJ = NAK(J)
         DO I = 1, NCORE
            IF (NP(I)/=NPJ .OR. NAK(I)/=NAKJ) CYCLE
            WRITE (ISTDE, *) 'lodcsh: The lists of core and', &
               ' peel subshells must form disjoint sets.'
            STOP
         END DO
      END DO

      WRITE (6, *) 'There are/is ', NW, ' relativistic subshells;'
!
!   Skip the header for the list of CSFs
!
      READ (NFILE, *)
!
!   To determine the number of electrons. This was done very much later
!   in the non-block mcp program( near the end of lodcsh). In the block
!   version, that will be used as a check to the value obtained below.
!
      READ (NFILE, '(A)', IOSTAT=IOS) STR
      CALL PRSRCN (STR, NCORE, IOCC, IERR)
      BACKSPACE (NFILE)                          ! return to the first CSF item
!             Number of electrons in the peel shells
      NELEC = SUM(IOCC(NCORE+1:NW))
!             Add the number of electrons in the core shells
      NELEC = NELEC + SUM(NKJ(:NCORE)+1)
      RETURN
      END SUBROUTINE LODCSH
