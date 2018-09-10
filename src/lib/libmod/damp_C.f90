      MODULE damp_C 
!     Arrays associated with damping solutions in the SCF process
!     ... ODAMP:  array for damping radial functions
!     ... CDAMP:  array for damping expansion coefficients
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNW
      REAL(DOUBLE), DIMENSION(NNNW) :: ODAMP 
      REAL(DOUBLE), DIMENSION(:), pointer :: cdamp
      END MODULE damp_C 
