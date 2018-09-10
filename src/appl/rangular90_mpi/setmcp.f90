!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
      SUBROUTINE SETMCP(NCORE, IDBLK, FILEHEAD) 
!                                                                      *
!   Open and check the  .mcp  files. File 30 stores the structure of   *
!   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
!   store V(0), V(1), ... .                                            *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN, OPENFL.                        *
!               [GENMCP]: GETINF.                                      *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 30 Jun 1998   *
!
!   Used by mcpvu, mcpmpivu
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mcp_C,     ONLY: KMAX
      USE orb_C,     ONLY: NW, NKJ
      USE hblock_C,  ONLY: NBLOCK, NCFBLK
      USE iounit_C,  ONLY: ISTDE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I 
      USE openfl_I 
      USE mpi_C,     ONLY: MYID, NPROCS, IERR
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     INTEGER, INTENT(IN) :: MYID 
!     INTEGER, INTENT(IN) :: NPROCS 
      INTEGER, INTENT(IN) :: NCORE 
      CHARACTER, INTENT(IN) :: FILEHEAD*(*) 
      CHARACTER, INTENT(IN) :: IDBLK(*)*8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, LNG, LCK, I 
      LOGICAL :: FOUND, FOUND1, GETYN, YES 
      CHARACTER :: CK*2 
!-----------------------------------------------
!   Determine KMAX; this is the number of  .mcp  files for the
!   two-electron integrals
 
      KMAX = 0 
      DO K = 1, NW 
         KMAX = MAX(KMAX,NKJ(K)) 
      END DO 
 
!   All files  mcp.xx  are UNFORMATTED;
 
      LNG = LEN_TRIM(FILEHEAD) 
      DO K = 30, 32 + KMAX 
         CALL CONVRT (K, CK, LCK) 
         CALL OPENFL (K, FILEHEAD(1:LNG)//'.'//CK(1:2), 'UNFORMATTED', &
            'UNKNOWN', IERR) 
         IF (IERR == 0) CYCLE  
         DO I = 30, K 
            CLOSE(I) 
         END DO 
         WRITE (ISTDE, *) 'Error when opening the mcp files' 
         STOP  
      END DO 
!
!   We want to know kmax before openning other mcp files (not mcp.30)
!   in rscf
!
      WRITE (30) NCORE, NBLOCK, KMAX 
      WRITE (30) (NCFBLK(I),I=1,NBLOCK) 
      WRITE (30) (IDBLK(I),I=1,NBLOCK) 
 
      DO K = 30, 32 + KMAX 
         WRITE (K) 'MCP', NBLOCK, MYID, NPROCS 
      END DO 
 
      RETURN  
      END SUBROUTINE SETMCP 
