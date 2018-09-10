!***********************************************************************
 
      SUBROUTINE SETMCP(NCORE, NBLKIN, IDBLK, FILEHEAD) 
!
!   Open, read, check the header of all .mcp files. Info for each
!   block is not accessed here.
!   File 30 stores the structure of H(DC) ;
!   file 31 stores the  T  coefficients;  files 32, 33, ...,
!   store V(0), V(1), ...
!
!   This version works in both serial/parallel(mpi) environment. The
!   difference comes from parameters filehead (in the arg list), myid
!   and nprocs (in common/mpi/).
!
!   The following items are read from mcp files. myid and nprocs
!   read from the files are checked against the ones in COMMON/mpi/
!
!       ncore, nblock, kmaxf             mcp.30 only
!       ncfblk(i), i=1, nblock)          mcp.30 only
!       idblk(i), i=1, nblock)             mcp.30 only
!
!       MCPLAB, nblock, myid, nprocs     all mcp files
!       nelec, ncf, nw               all mcp files
!       DIAG, ICCUT, LFORDR                all mcp files
!
!   Call(s) to: [LIB92]: CONVRT, OPENFL.
!
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!   Modified by Xinghong He               Last revision: 06 Aug 1998   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:27:31   1/ 5/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE memory_man
      USE DEF_C
      USE FOPARM_C
      USE MCPA_C 
      USE MCPB_C
      USE mpi_C
      USE orb_C
      USE iounit_C
      USE hblock_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I 
      USE openfl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE 
      INTEGER , INTENT(IN) :: NBLKIN 
      CHARACTER , INTENT(IN) :: FILEHEAD*(*) 
      CHARACTER  :: IDBLK(*)*8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENFH, IERROR, IOS, I, K, LCK, MYIDD, NPROCSS 
      LOGICAL :: FOUND, FOUND1 
      CHARACTER :: CK*2, MCPLAB*3 
      CHARACTER(LEN=120) :: FILNAM
!-----------------------------------------------
!
      LENFH = LEN_TRIM(FILEHEAD) 
 
      FILNAM = FILEHEAD(1:LENFH)//'.30' 
      OPEN(30, FILE=FILNAM, FORM='UNFORMATTED', STATUS='OLD', IOSTAT=IERROR, &
         POSITION='asis') 
 
!  Parameter ierror carries through mcp.30,... mcp.kmax
 
      READ (30, IOSTAT=IOS) NCORE, NBLOCK, KMAXF 
      IERROR = IERROR + ABS(IOS) 
      IF (NBLOCK>NBLKIN .OR. NBLOCK<1) THEN 
         WRITE (ISTDE, *) 'setmcp: nblock = ', NBLOCK 
         STOP  
      ENDIF 
 
!cjb allocate ncfblk(0:*)
!cjb  CALL ALLOC (NCFBLK, NBLOCK + 1, 'NCFBLK', 'SETMCP')
      CALL ALLOC (NCFBLK, 0, NBLOCK , 'NCFBLK', 'SETMCP')
!cjb
      NCFBLK(0) = 0 
 
      READ (30, IOSTAT=IOS) (NCFBLK(I),I=1,NBLOCK) 
      IERROR = IERROR + ABS(IOS) 
      READ (30, IOSTAT=IOS) (IDBLK(I),I=1,NBLOCK) 
      IERROR = IERROR + ABS(IOS) 
 
!  Look for other mcp files
 
      FOUND = .TRUE. 
      DO K = 31, 32 + KMAXF 
         CALL CONVRT (K, CK, LCK) 
         FILNAM = FILEHEAD(1:LENFH)//'.'//CK(1:2) 
         INQUIRE(FILE=FILNAM, EXIST=FOUND1) 
         FOUND = FOUND .AND. FOUND1 
      END DO 
 
      IF (.NOT.FOUND) THEN 
         WRITE (ISTDE, *) 'The mcp files do not exist' 
         STOP  
      ENDIF 
 
!  Open the files; check file headers
 
      DO K = 30, 32 + KMAXF 
 
         IF (K /= 30) THEN 
            CALL CONVRT (K, CK, LCK) 
            FILNAM = FILEHEAD(1:LENFH)//'.'//CK(1:2) 
            CALL OPENFL (K, FILNAM, 'UNFORMATTED', 'OLD', IERROR) 
         ENDIF 
 
         READ (K, IOSTAT=IOS) MCPLAB, NBLOCK, MYIDD, NPROCSS 
 
         IERROR = IERROR + ABS(IOS) 
         IF (MYID/=MYIDD .OR. NPROCS/=NPROCSS) THEN 
            WRITE (ISTDE, *) 'mcp files were generated under different', &
               ' processor configuration.' 
            STOP  
         ENDIF 
 
         IF (MCPLAB /= 'MCP') THEN 
            WRITE (ISTDE, *) 'Not a sorted GRASP92 MCP File;' 
            IERROR = IERROR + 1 
         ENDIF 
 
         READ (K, IOSTAT=IOS) NELEC, NCF, NW 
         IERROR = IERROR + ABS(IOS) 
         READ (K, IOSTAT=IOS) DIAG, ICCUT, LFORDR 
         IERROR = IERROR + ABS(IOS) 
 
         IF (IERROR == 0) CYCLE  
         WRITE (ISTDE, *) 'setmcp: Error accumulated , stopping...' 
         DO I = 30, K 
            CLOSE(I) 
         END DO 
         STOP  
      END DO 
 
      RETURN  
      END SUBROUTINE SETMCP 
