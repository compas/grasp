!***********************************************************************
!                                                                      *
      SUBROUTINE EVCOUT
!                                                                      *
!   Routine for printing eigenvectors.                                 *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man
      USE eigv_C
      USE orb_C,           ONLY: ncf, nw, iqa
      USE prnt_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V A R I A B L E S
!-----------------------------------------------
      INTEGER :: I, J1, J2, nvex, jp, l, m
      INTEGER, DIMENSION(:), pointer :: ILOC
!-----------------------------------------------
!
!   Header
!
      WRITE (24,300)
!
!   Determine which eigenvectors are to be printed out
!
      NVEX = NVEC
      CALL ALLOC (ILOC,NVEC,'ILOC', 'EVCOUT')
      DO 1 I = 1,NVEC
         ILOC(I) = IVEC(I)
    1 CONTINUE
!
!   There are eight columns across the width of a page
!
      JP = 8
!
!   Set up for the first set of columns
!
      J1 = 1
      J2 = JP
!
!   Loop over sets of columns
!
    3 IF (J2 .GT. NVEX) J2 = NVEX
!
!   Write out the column numbers; skip a line
!
      WRITE (24,301) (ILOC(L),L = J1,J2)
      WRITE (24,302)
!
!   Write out the rows
!
      DO 4 M = 1,NCF
!
         WRITE (24,303) (EVEC(M+(L-1)*NCF),L = J1,J2)
!
!   Skip a line after every tenth row
!
         IF (MOD(M,10) .EQ. 0) WRITE (24,302)
!
    4 CONTINUE
!
!   Set up for the next set of columns
!
      WRITE (24,303)
      IF (J2 .LT. NVEX) THEN
         J1 = J1+JP
         J2 = J2+JP
         GOTO 3
      ENDIF
!
!   Deallocate storage for ILOC
!
      CALL DALLOC (ILOC, 'ILOC', 'EVCOUT')
!
      RETURN
!
  300 FORMAT (//' Eigenvectors:'/)
  301 FORMAT (1X,8(I8,7X))
  302 FORMAT (1X)
  303 FORMAT (1X,1P,8D15.7)
!
      END
