!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      SUBROUTINE SETMCPmpi(NCORE, IDBLK, FILEHEAD) 
!                                                                      *
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial *
! and MPI programs.                                                    *
!                                                                      *
!   Written by Xinghong He                Last revision: 30 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!cjb mpi_C
!     USE mpi_C,    ONLY: MYID, NPROCS, ierr
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE iccu_C
      USE mcp_C
      USE orb_C
      USE def_C
      USE hblock_C
!cjb mpi_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setmcp_I 
!     USE setmcpmpi_I
      USE getinf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      INTEGER , INTENT(IN) :: MYID
!      INTEGER , INTENT(IN) :: NPROCS 
       INTEGER , INTENT(IN) :: NCORE 
      CHARACTER , INTENT(IN) :: FILEHEAD*(*) 
       CHARACTER , INTENT(IN) :: IDBLK(*)*8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!cjb mpi_C
!cjb  INTEGER :: K, ierr
      INTEGER :: K
!-----------------------------------------------

      CALL SETMCP (ncore, idblk, filehead)

! DIAG, ICCUT, LFORDR are set in GETINF

      IF (myid .EQ. 0) THEN
         CALL GETINF
      ENDIF

      CALL MPI_Bcast (ICCUT(1),100,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (DIAG,  1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (LFORDR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      DO K = 30, 32+KMAX
         WRITE (K) NELEC,NCF,NW
         WRITE (K) DIAG,ICCUT(1),LFORDR
      ENDDO

      DO K = 30, 32+KMAX
         CLOSE (K)
      ENDDO
 
      RETURN  
      END SUBROUTINE SETMCPmpi
