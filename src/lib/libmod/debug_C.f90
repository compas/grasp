!
!***********************************************************************
!                                                                      *
      MODULE debug_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      INTEGER :: IBUG1, IBUG2, IBUG3, IBUG4, IBUG5, IBUG6 
      LOGICAL, DIMENSION(5) :: LDBPA 
      LOGICAL, DIMENSION(5) :: LDBPG 
      LOGICAL, DIMENSION(30) :: LDBPR 
      REAL(DOUBLE) :: cutoff  ! used by bioscl
      END MODULE debug_C 
