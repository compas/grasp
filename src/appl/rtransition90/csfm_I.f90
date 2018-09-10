      MODULE csfm_I   
      INTERFACE
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
      SUBROUTINE CSFM (ASFA,ASFB,LEV1,LEV2)
        USE vast_kind_param, ONLY:  DOUBLE    
        REAL(DOUBLE), INTENT(OUT) :: asfa, asfb
        INTEGER, INTENT(IN) :: lev1, lev2
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
