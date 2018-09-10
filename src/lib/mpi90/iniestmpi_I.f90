      MODULE INIESTmpi_I
      INTERFACE
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE INIESTmpi (NMAX, NCF, NIV, BASIS, HMX, JCOL, IROW)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER, INTENT(IN) :: NMAX
      INTEGER, INTENT(IN) :: NCF
      INTEGER  :: NIV
      INTEGER, INTENT(IN) :: JCOL(0:*)
      INTEGER, INTENT(IN) :: IROW(*)
      real(double)  :: BASIS(*)
      real(double), INTENT(IN) :: HMX(*)
      END SUBROUTINE
      END INTERFACE
      END MODULE
