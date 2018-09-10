!***********************************************************************
!                                                                      *
      SUBROUTINE RADFILE(NAME) 
!                                                                      *
!  This subroutine outputs the transformed radial orbitals             *
!                                                                      *
!  Written by Per Jonsson                                              *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE biorb_C
      USE grid_C
      USE orb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER, INTENT(IN) :: NAME(2)*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, I 
!-----------------------------------------------
!
      J = INDEX(NAME(1),' ') 
      OPEN(UNIT=30, FILE=NAME(1)(1:J-1)//'.bw', FORM='UNFORMATTED', STATUS=&
         'UNKNOWN', POSITION='asis') 
 
      WRITE (30) 'G92RWF' 
      WRITE (*, *) 'NWII', NWII 
      DO K = 1, NWII 
         WRITE (30) NPII(K), NAKII(K), EII(K), MFII(K) 
         WRITE (30) PZII(K), (PFII(I,K),I=1,MFII(K)), (QFII(I,K),I=1,MFII(K)) 
         WRITE (30) (R(I),I=1,MFII(K)) 
      END DO 
 
      CLOSE(30) 
 
      J = INDEX(NAME(2),' ') 
      OPEN(UNIT=30, FILE=NAME(2)(1:J-1)//'.bw', FORM='UNFORMATTED', STATUS=&
         'UNKNOWN', POSITION='asis') 
 
      WRITE (30) 'G92RWF' 
      WRITE (*, *) 'NWFF', NWFF 
      DO K = 1, NWFF 
         WRITE (30) NPFF(K), NAKFF(K), EFF(K), MFFF(K) 
         WRITE (30) PZFF(K),(PFFF(I,K),I=1,MFFF(K)),(QFFF(I,K),I=1,MFFF(K)) 
         WRITE (30) (R(I),I=1,MFFF(K)) 
      END DO 
 
      CLOSE(30) 
 
      RETURN  
      END SUBROUTINE RADFILE 
