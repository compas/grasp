!***********************************************************************
!                                                                      *
      SUBROUTINE SETSDA(OUTSDAdummy,NNONZ,LPRINT,NB,MYID,NPROCS,FHEAD) 
!                                                                      *
!   This routine examines lists                                        *
!                               (IC,IR,npos)                           *
!   to set up the array  IENDC  required by the Davidson eigensolver   *
!   of Stathopoulos and Fischer.                                       *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
!   Modified for block interactions by C. F. Fischer        May 1997   *
!   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
!
!  Currently shared by mcpblk, mcpmpi
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:52:18   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE orb_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE outsdampi_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      EXTERNAL OUTSDAdummy
      INTEGER, INTENT(IN) :: NNONZ
      INTEGER , INTENT(IN) :: NB 
      INTEGER , INTENT(IN) :: MYID 
      INTEGER , INTENT(IN) :: NPROCS 
      LOGICAL  :: LPRINT 
      CHARACTER , INTENT(IN) :: FHEAD*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MB, IEND, ICLAST, I, IC, NPOS, IERR 
      INTEGER, DIMENSION(:), pointer :: IENDC
      INTEGER, DIMENSION(:), pointer :: IROW 
      CHARACTER :: MCPLAB*3 
!-----------------------------------------------
 
      IF (MYID == 0) WRITE (6, *) &
         'Analysing sparse matrix array definition file ...', 30 
 
      READ (30) MCPLAB, MB 
      IF (NB /= MB) THEN 
         WRITE (ISTDE, *) 'setsda: nb = ', NB, '.NE. mb (=', MB, ')' 
         STOP  
      ENDIF 
!
!cjb           PRINT *, ' setsda.f90 FHEAD '
!cjb           PRINT *, ' setsda.f90 FHEAD = ', FHEAD
!   Allocate storage for IENDC(0:NCF)
!
      CALL ALLOC (IENDC, 0,  NCF,'IENDC','SETSDA'  ) 
      CALL ALLOC (IROW, NNONZ, 'IROW', 'SETSDA' ) 
!
!   Analyse data on file 30; set up IENDC and IROW
! In multiprocessor environment, iendc of each node will have the
! same length (ncf+1); but will have its own part filled. irow is
! local, and its length is determined by the local parameter nnonz.
 
      IEND = 0 
      ICLAST = 0 
      DO I = 1, NNONZ 
         READ (30) IC, IROW(I), NPOS 
         IF (IC /= ICLAST) THEN 
            IENDC(ICLAST) = IEND 
            ICLAST = IC 
         ENDIF 
         IEND = NPOS 
      END DO 
!xhh - changed to suits MPI environment as well
!      IENDC(NCF) = IEND
      IENDC(IC) = IEND 
!
!   Sorting complete; rewrite to mcpXXX.30 file
!
!cjb           PRINT *, ' before FHEAD '
!cjb           PRINT *, ' FHEAD = ', FHEAD
!cjb           PRINT *, ' before OPEN '
      OPEN(29, FILE=FHEAD//'.30', STATUS='OLD', FORM='UNFORMATTED', IOSTAT=IERR&
         , POSITION='APPEND') 
!cjb           PRINT *, ' after OPEN '
      IF (IERR /= 0) THEN 
         WRITE (ISTDE, *) ' Error when opening the file mcp.30' 
         STOP  
      ENDIF 
 
      WRITE (29) 'MCP', NB, NCF 
      WRITE (29) NNONZ 
      WRITE (29) (IENDC(I),I=MYID + 1,NCF,NPROCS), (IROW(I),I=1,NNONZ) 
      CLOSE(29) 
!cjb
!
!   Deallocate storage
!
      CALL DALLOC (IENDC, 'IENDC', 'SETSDA') 
      CALL DALLOC (IROW, 'IROW', 'SETSDA') 
!
      CALL OUTSDAdummy (LPRINT, NNONZ, NCF)
      RETURN  
      END SUBROUTINE SETSDA 
