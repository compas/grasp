!************************************************************************
!*									*
      SUBROUTINE GETMIXBLOCK(NAME, NCI)
!*									*
!*    Reads mixing coefficient file from block-structured format	*
!*									*
!* Note:								*
!*     eav is not compatible with the non-block version if some blocks  *
!*     were not diagonalized						*
!*									*
!*     This is a modified version of cvtmix.f				*
!*     written by Per Jonsson, September 2003				*
!*									*
!*   Translated by Pacific-Sierra Research 77to90  4.3E  18:32:57 1/6/07*  
!*   Modified by Charlotte Froese Fischer 				*
!*                     Gediminas Gaigalas  11/01/17			*
!*   Modified by Wenxian Li F77 to F90 12/28/18                         * 
!************************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE memory_man
      USE def_C
      USE EIGV_C 
      USE orb_C
      USE prnt_C
      USE syma_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)           :: NCI 
      CHARACTER(LEN=24), INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IERR, IOS, NCFTOT, NVECTOT, NVECSIZ, NBLOCK, I, NVECPAT, &
         NCFPAT, NVECSIZPAT, NEAVSUM, JB, NB, NCFBLK, NEVBLK, IATJP, IASPA, J 
      REAL(DOUBLE) :: EAVSUM 
      CHARACTER :: FILNAM*256, FORM*11, G92MIX*6, STATUS*3 
!-----------------------------------------------
!
!   The  .mix  file is UNFORMATTED; it must exist
!
      K = INDEX(NAME,' ') 
      IF (NCI == 0) THEN 
         FILNAM = NAME(1:K-1)//'.cm' 
      ELSE 
         FILNAM = NAME(1:K-1)//'.m' 
      ENDIF 
      FORM = 'UNFORMATTED' 
      STATUS = 'OLD' 
!
      CALL OPENFL (25, FILNAM, FORM, STATUS, IERR) 
      IF (IERR == 1) THEN 
         WRITE (ISTDE, *) 'Error when opening', FILNAM 
         STOP  
      ENDIF 
!
!   Check the header of the file; if not as expected, try again
!
      READ (25, IOSTAT=IOS) G92MIX 
      IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN 
         WRITE (ISTDE, *) 'Not a GRASP92 MIXing Coefficients File;' 
         CLOSE(25) 
         STOP  
      ENDIF 
 
      READ (25) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK 
      WRITE (*, *) '   nelec  = ', NELEC 
      WRITE (*, *) '   ncftot = ', NCFTOT 
      WRITE (*, *) '   nw     = ', NW 
      WRITE (*, *) '   nblock = ', NBLOCK 
      WRITE (*, *) 
 
!***********************************************************************
! Allocate memory for old format data
!***********************************************************************
 
      CALL ALLOC (EVAL, NVECTOT, 'EVAL', 'GETMIXBLOCK') 
      CALL ALLOC (EVEC, NCFTOT*NVECTOT, 'EVEC', 'GETMIXBLOCK') 
      CALL ALLOC (IVEC, NVECTOT, 'IVEC', 'GETMIXBLOCK') 
      CALL ALLOC (IATJPO, NVECTOT, 'IATJPO', 'GETMIXBLOCK') 
      CALL ALLOC (IASPAR, NVECTOT, 'IASPAR', 'GETMIXBLOCK') 
 
!***********************************************************************
! Initialize mixing coefficients to zero; others are fine
!***********************************************************************
      EVEC(:NVECTOT*NCFTOT) = 0.D0 
 
!***********************************************************************
! Initialize counters and sum registers
!
!    nvecpat:    total number of eigenstates of the previous blocks
!    ncfpat:     total number of CSF of the previous blocks
!    nvecsizpat: vector size of the previous blocks
!    eavsum:     sum of diagonal elements of the previous blocks where
!                at least one eigenstate is calculated
!    neavsum:    total number CSF contributing to eavsum
!***********************************************************************
 
      NVECPAT = 0 
      NCFPAT = 0 
      NVECSIZPAT = 0 
      NEAVSUM = 0 
      EAVSUM = 0.D0 
 
      WRITE (*, *) '  block     ncf     nev    2j+1  parity' 
      DO JB = 1, NBLOCK 
 
         READ (25) NB, NCFBLK, NEVBLK, IATJP, IASPA 
         WRITE (*, '(5I8)') NB, NCFBLK, NEVBLK, IATJP, IASPA 
         IF (JB /= NB) STOP 'jb .NE. nb' 
 
         IF (NEVBLK > 0) THEN 
 
            READ (25) (IVEC(NVECPAT + I),I=1,NEVBLK) 
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
            IATJPO(NVECPAT+1:NEVBLK+NVECPAT) = IATJP 
            IASPAR(NVECPAT+1:NEVBLK+NVECPAT) = IASPA 
 
            READ (25) EAV, (EVAL(NVECPAT+I),I=1,NEVBLK) 
 
!           ...Construct the true energy by adding up the average
            EVAL(NVECPAT+1:NEVBLK+NVECPAT) = EVAL(NVECPAT+1:NEVBLK+NVECPAT) + &
               EAV 
!           ...For overal (all blocks) average energy
            EAVSUM = EAVSUM + EAV*NCFBLK 
            NEAVSUM = NEAVSUM + NCFBLK 
 
            READ (25) ((EVEC(NVECSIZPAT+NCFPAT+I+(J-1)*NCFTOT),I=1,NCFBLK),J=1,&
               NEVBLK) 
         ENDIF 
 
         NVECPAT = NVECPAT + NEVBLK 
         NCFPAT = NCFPAT + NCFBLK 
         NVECSIZPAT = NVECSIZPAT + NEVBLK*NCFTOT 
 
      END DO 
 
!     ...Here eav is the average energy of the blocks where at least
!        one eigenstate is calculated. It is not the averge of the
!        total Hamiltonian.
 
      EAV = EAVSUM/NEAVSUM 
 
      IF (NCFTOT /= NEAVSUM) WRITE (6, *) &
         'Not all blocks are diagonalized --- Average E ', 'not correct'
 
!     ...Substrct the overal average energy
      EVAL(:NVECTOT) = EVAL(:NVECTOT) - EAV 
 
      CLOSE(25) 
 
      NCF = NCFTOT 
      NVEC = NVECTOT 
 
      RETURN  
      END SUBROUTINE GETMIXBLOCK 
