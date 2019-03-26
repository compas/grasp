!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBGmpi(DBGFILE)
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 10 Dec 1992   *
!   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:22:29   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE DEBUG_C
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) , INTENT(IN) :: DBGFILE
!-----------------------------------------------
!
      IF (myid .EQ. 0) THEN
         CALL setdbg (dbgfile)
      ENDIF

      CALL MPI_Bcast (ldbpa, 5, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ldbpg, 5, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ldbpr,30, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      RETURN
      END SUBROUTINE SETDBGmpi
