      MODULE mcp_I
      INTERFACE
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE mcp (nb, RESTRT, myid, nprocs, fhead)
      INTEGER  :: NB
      INTEGER  :: MYID
      INTEGER  :: NPROCS
      LOGICAL, INTENT(IN) :: RESTRT
      CHARACTER(len=*), INTENT(IN) :: FHEAD
      END SUBROUTINE
      END INTERFACE
      END MODULE
