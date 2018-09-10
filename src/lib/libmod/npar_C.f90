!
!***********************************************************************
!                                                                      *
      MODULE npar_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: NPARM 
      REAL(DOUBLE), DIMENSION(2) :: PARM 
      END MODULE npar_C 

      MODULE nsmdat_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:36:34  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      REAL(DOUBLE) :: SQN,DMOMNM,QMOMB
      REAL(DOUBLE) :: HFSI, HFSD, HFSQ 
      REAL(DOUBLE) :: SMSI, SMSD, SMSQ 
      END MODULE nsmdat_C 
