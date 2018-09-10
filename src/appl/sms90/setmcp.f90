!***********************************************************************
!                                                                      *
      SUBROUTINE SETMCP(AVAIL) 
!                                                                      *
!   Open and check the  .mcp  files. File 30 stores the structure of   *
!   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
!   store V(0), V(1), ... .                                            *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, LENGTH, OPENFL.                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE def_C
      USE iccu_C
      USE mcp_C,   KMAXF=>KMAX
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I 
      USE openfl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(OUT) :: AVAIL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, LCK, LFN, IERR, IOS, NELECT, NCFT, NWT, ICCUTT, I 
      LOGICAL :: DIAGT, FOUND, FOUND1, LFORDT 
      CHARACTER :: FILNAM*256, FULNAM*256, DEFNAM*3, FORM*11, SRTLAB*8, MCPLAB&
         *3, STATUS*3, CK*2 
!-----------------------------------------------
!
!   Determine KMAXF; this is one less than the number of  .mcp
!   files for the two-electron integrals
!
      AVAIL = .TRUE. 
      KMAXF = 0 
      DO K = 1, NW 
         KMAXF = MAX(KMAXF,NKJ(K)) 
      END DO 
!
!   All files  grasp92.mcp.xx  are UNFORMATTED; they must exist
!
      FORM = 'UNFORMATTED' 
      DEFNAM = 'mcp' 
      STATUS = 'OLD' 
!
!   Look for  grasp92.mcp.30 , ...
!
      FOUND = .TRUE. 
      DO K = 30, 32 + KMAXF 
         CALL CONVRT (K, CK, LCK) 
         INQUIRE(FILE=DEFNAM//'.'//CK(1:2), EXIST=FOUND1) 
         FOUND = FOUND .AND. FOUND1 
      END DO 
!
      IF (FOUND) THEN 
         FILNAM = DEFNAM 
      ELSE 
         WRITE (6, *) 'The mcp files does not exist' 
         AVAIL = .FALSE. 
         RETURN  
      ENDIF 
!
!   Open the files; check file headers
!
      LFN = LEN_TRIM(FILNAM) 
      DO K = 30, 32 + KMAXF 
         CALL CONVRT (K, CK, LCK) 
         FULNAM = FILNAM(1:LFN)//'.'//CK(1:2) 
         CALL OPENFL (K, FULNAM, FORM, STATUS, IERR) 
         IF (IERR == 0) THEN 
            READ (K, IOSTAT=IOS) MCPLAB, SRTLAB 
            IF (IOS/=0 .OR. MCPLAB/='MCP' .OR. SRTLAB/='  SORTED') THEN 
               WRITE (6, *) 'Not a sorted GRASP92 MCP File;' 
               IERR = IERR + 1 
            ENDIF 
         ENDIF 
         IF (IERR == 0) THEN 
            READ (K) NELECT, NCFT, NWT 
            IF (NELECT/=NELEC .OR. NCFT/=NCF .OR. NWT/=NW) THEN 
               WRITE (6, *) 'Sorted GRASP92 MCP File not appropriate' 
               WRITE (6, *) '  to Configuration Symmetry List;' 
               IERR = IERR + 1 
            ENDIF 
            IF (K == 30) THEN 
               READ (K) DIAG, ICCUT(1), LFORDR 
            ELSE 
               READ (K) DIAGT, ICCUTT, LFORDT 
               IF ((DIAGT .NEQV. DIAG) .OR. ICCUTT/=ICCUT(1) .OR. (LFORDT .NEQV. &
                  LFORDR)) THEN 
                  WRITE (6, *) 'Sorted GRASP92 MCP Files are not' 
                  WRITE (6, *) ' consistent;' 
                  IERR = IERR + 1 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (IERR == 0) CYCLE  
         DO I = 30, K 
            CLOSE(I) 
         END DO 
         AVAIL = .FALSE. 
         RETURN  
      END DO 
!
      RETURN  
      END SUBROUTINE SETMCP 
