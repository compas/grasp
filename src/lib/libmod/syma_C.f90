!
!***********************************************************************
!                                                                      *
      MODULE syma_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
      INTEGER, PARAMETER::nblk0 = 50    ! Maximum number of blocks
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER, DIMENSION(:), pointer :: iatjpo, iaspar
      INTEGER, DIMENSION(:), pointer :: iatjpoff, iasparff
      INTEGER, DIMENSION(:), pointer :: iatjpoii, iasparii
      INTEGER, DIMENSION(nblk0) :: jpgg
      END MODULE syma_C 
