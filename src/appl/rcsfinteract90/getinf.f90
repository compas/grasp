!***********************************************************************
!                                                                      *
      SUBROUTINE GETINF
!                                                                      *
!   Interactively determines data governing the generation of MCP co-  *
!   efficients.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!   Modified Xinghong He                  Last revision: 30 Jun 1998   *
!                                                                      *
!   File shared (hard link) by mcpmpi, mcpblk
!
!   Updated to treat ICCUT for block
!
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE CONS_C
      USE orb_C
      USE iccu_C
      USE default_C
      USE iounit_C
      USE mcp_C
      USE hblock_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------

!-----------------------------------------------------------------------
!
!   Determine the physical effects specifications
!
!     Commenting out the EAL option
!     IF (NDEF .NE. 0) THEN
!        WRITE  (istde,*) 'Generate MCP coefficients only for'
!    & , ' diagonal matrix elements? '
!        WRITE (istde,*) '(This is appropriate to (E)AL calculation):'
!        DIAG = GETYN ()
!     ELSE
!        DIAG = .FALSE.
!     ENDIF
      DIAG = .FALSE.
      IF (DIAG) THEN
         LFORDR = .FALSE.
         do i = 1,100
            ICCUT(i) = 0
         end do
      ELSE
         IF (NDEF .NE. 0) THEN
!            WRITE (istde,*) 'Treat contributions of some CSFs', &
!                            ' as first-order perturbations?'
!            LFORDR = GETYN ()
            LFORDR = .TRUE.
         ELSE
            LFORDR = .FALSE.
         ENDIF
         IF (LFORDR) THEN
            WRITE (istde,*) 'The contribution of CSFs 1 -- ICCUT will',&
                            ' be treated variationally;'
            WRITE (istde,*) 'the remainder perturbatively; enter ICCUT:'
            do i = 1,nblock
               write(istde,*) 'Give ICCUT for block',i
               READ *, ICCUT(i)
               write(739,*) ICCUT(i), '! ICCUT FOR BLOCK',i
!    1          READ *, ICCUT(i)
!               IF ((ICCUT(i).LE.1).OR.(ICCUT(i).GE.ncfblk(i))) THEN
!                  WRITE (istde,*) 'GETINF: ICCUT must be greater than 1', &
!                                  ' and less than ',ncfblk(i)
!                  WRITE (istde,*) ' please reenter ICCUT:'
!                  GOTO 1
!               ENDIF
            end do
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE GETINF
