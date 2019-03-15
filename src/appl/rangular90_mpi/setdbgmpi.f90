!***********************************************************************
!                                                                      *
      SUBROUTINE SETDBGmpi (DEBUG, fullname)
!                                                                      *
!   This routine calls setdbg to open the  .dbg  file and sets the     *
!   arrays that control debug printout from the GENMCP program.        *
!                                                                      *
!   Written by Xinghong He                  Last update: 29 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:11:16  12/23/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE MPI_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setdbg_I
      USE DEBUG_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(out) :: debug
      Character(LEN=*) :: fullname
!-----------------------------------------------
!
     IF (myid == 0) THEN
         CALL setdbg (debug, fullname)
      ENDIF
      CALL MPI_Bcast (debug,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ldbpa,5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ldbpg,5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE SETDBGmpi
