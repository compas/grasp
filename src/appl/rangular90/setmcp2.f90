!***********************************************************************
      SUBROUTINE SETMCP2(MYID, NPROCS, NCORE, IDBLK, FILEHEAD)
!
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
! and MPI programs.
!
!   Written by Xinghong He                Last revision: 30 Jun 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE def_C
      USE iccu_C
      USE mcp_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setmcp_I
      USE getinf_I
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      INTEGER  :: NCORE
      CHARACTER  :: FILEHEAD*(*)
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K

      CALL SETMCP (MYID, NPROCS, NCORE, IDBLK, FILEHEAD)

! DIAG, ICCUTblk, LFORDR are set in GETINF

      IF (MYID == 0) CALL GETINF

      DO K = 30, 32 + KMAX
         WRITE (K) NELEC, NCF, NW
         WRITE (K) DIAG, ICCUT(1), LFORDR
      END DO

      DO K = 30, 32 + KMAX
         CLOSE(K)
      END DO

      RETURN
      END SUBROUTINE SETMCP2
