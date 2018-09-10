!***********************************************************************
!                                                                      *
      SUBROUTINE TESTMIX 
!                                                                      *
!   This routine checks the mixing coefficients                        *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE biorb_C
      USE def_C
      USE eigv_C
      USE prnt_C
      USE syma_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J 
!-----------------------------------------------
!
!
 
      WRITE (*, *) ' ****************' 
      WRITE (*, *) ' Entering testmix' 
      WRITE (*, *) ' ****************' 
      WRITE (*, *) 
      WRITE (*, *) 'Initial state' 
      WRITE (*, *) 'EVALII', (EAVII + EVALII(I),I=1,NVECII) 
      WRITE (*, *) 'NELECII,NCFII,NWII,NVECMXII', NELECII, NCFII, NWII, NVECMXII
      WRITE (*, *) NVECII 
      WRITE (*, *) (IVECII(I),I=1,NVECII) 
      WRITE (*, *) (IATJPOII(I),IASPARII(I),I=1,NVECII) 
      WRITE (*, *) ((EVECII(I + (J - 1)*NCFII),I=1,NCFII),J=1,NVECII) 
 
      WRITE (*, *) 'Final state' 
      WRITE (*, *) 'EVALFF', (EAVFF + EVALFF(I),I=1,NVECFF) 
      WRITE (*, *) 'NELECFF,NCFFF,NWFF,NVECMXFF', NELECFF, NCFFF, NWFF, NVECMXFF
      WRITE (*, *) NVECFF 
      WRITE (*, *) (IVECFF(I),I=1,NVECFF) 
      WRITE (*, *) (IATJPOFF(I),IASPARFF(I),I=1,NVECFF) 
      WRITE (*, *) ((EVECFF(I + (J - 1)*NCFFF),I=1,NCFFF),J=1,NVECFF) 
      WRITE (*, *) 
      WRITE (*, *) ' ***************' 
      WRITE (*, *) ' Leaving testmix' 
      WRITE (*, *) ' ***************' 
 
      RETURN  
      END SUBROUTINE TESTMIX 
