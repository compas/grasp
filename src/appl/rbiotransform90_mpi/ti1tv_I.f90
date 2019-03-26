      MODULE ti1tv_I
      INTERFACE
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      SUBROUTINE TI1TV(CIIN,NCSF,NCIV,I,L,T,NSHL,CIOUT,NTESTG)
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ncsf, nciv, nshl, i, L, ntestg
      REAL(DOUBLE), DIMENSION(ncsf, nciv), INTENT(IN) :: ciin
      REAL(DOUBLE), DIMENSION(nshl), INTENT(IN) :: t
      REAL(DOUBLE), DIMENSION(ncsf, nciv), INTENT(OUT) :: ciout
      END SUBROUTINE
      END INTERFACE
      END MODULE
