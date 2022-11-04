!***********************************************************************
!                                                                      *
      SUBROUTINE READMIX(NAME, INPCI, INIT)
!                                                                      *
!   Open and read the mixing coefficent files                          *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:38   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE biorb_C
      USE def_C
      USE orb_C
      USE eigv_C
      USE prnt_C
      USE syma_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: INPCI
      INTEGER , INTENT(IN) :: INIT
      CHARACTER  :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IB, IATJP, IASPA, I, J
      CHARACTER :: G92MIX*6
!-----------------------------------------------
!
      IF (INIT == 1) THEN
!
!   Read the initial state mixing file
!
         READ (68) IB, NCFII, NVECII, IATJP, IASPA
         CALL ALLOC (EVALII, NVECII, 'EVALII', 'READMIX')
         CALL ALLOC (EVECII, NCFII*NVECII, 'EVECII', 'READMIX')
         CALL ALLOC (IVECII, NVECII,'IVECII', 'READMIX' )
         CALL ALLOC (IATJPOII, NVECII, 'IATJPOII', 'READMIX')
         CALL ALLOC (IASPARII, NVECII, 'ISPARII', 'READMIX')
         READ (68) (IVECII(I),I=1,NVECII)
         IATJPOII(:NVECII) = IATJP
         IASPARII(:NVECII) = IASPA
         READ (68) EAVII, (EVALII(I),I=1,NVECII)

         READ (68) ((EVECII(I + (J - 1)*NCFII),I=1,NCFII),J=1,NVECII)
!
!     CLOSE(68)
!
      ELSE
!
!   Read the final state mixing file
!

         READ (78) IB, NCFFF, NVECFF, IATJP, IASPA
         CALL ALLOC (EVALFF, NVECFF, 'EVALFF', 'READMIX')
         CALL ALLOC (EVECFF, NCFFF*NVECFF, 'EVECFF', 'READMIX')
         CALL ALLOC (IVECFF, NVECFF,'IVECFF', 'READMIX' )
         CALL ALLOC (IATJPOFF, NVECFF, 'IATJPOFF', 'READMIX')
         CALL ALLOC (IASPARFF, NVECFF, 'ISPARFF', 'READMIX')
         READ (78) (IVECFF(I),I=1,NVECFF)
         IATJPOFF(:NVECFF) = IATJP
         IASPARFF(:NVECFF) = IASPA
         READ (78) EAVFF, (EVALFF(I),I=1,NVECFF)

         READ (78) ((EVECFF(I + (J - 1)*NCFFF),I=1,NCFFF),J=1,NVECFF)
!
!   Close the initial state mixing  file
!
!      CLOSE(78)
      ENDIF

      RETURN
      END SUBROUTINE READMIX
