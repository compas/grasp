      MODULE cslhmpi_I
      INTERFACE
!
      SUBROUTINE cslhmpi (NAME, NCORE, NBLKIN, IDBLK)
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      CHARACTER (LEN = *) :: NAME
      INTEGER :: NCORE
      INTEGER, INTENT(IN) :: NBLKIN
      CHARACTER (LEN = 8), DIMENSION(*) :: IDBLK
      END SUBROUTINE
      END INTERFACE
      END MODULE
