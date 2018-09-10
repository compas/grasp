!
!***********************************************************************
!                                                                      *
      MODULE prnt_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NVEC, NVECMX 
      INTEGER :: NVECFF, NVECMXFF 
      INTEGER :: NVECII, NVECMXII 
      REAL(DOUBLE) :: PNIVECII 
      INTEGER, DIMENSION(:), pointer :: ivec, ivecff, ivecii
      END MODULE prnt_C 
