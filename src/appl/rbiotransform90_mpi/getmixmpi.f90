!***********************************************************************
!                                                                      *
      SUBROUTINE GETMIX(NAME, INPCI, IBLK)
!                                                                      *
!   Open, check, load data from and close the rscf.mix file.           *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, LENGTH, OPENFL.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:21:54   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE def_C
      USE orb_C, ONLY: ncf, nw, iqa
      USE EIGV_C
      USE PRNT_C
      USE SYMA_C
      USE BLK_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: INPCI
      INTEGER, INTENT(IN) :: IBLK
      CHARACTER, INTENT(IN) :: NAME*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IK, J, IOS, NB, IATJP, IASPA, I
      CHARACTER :: G92MIX*6
!-----------------------------------------------

      IK = 30
      IF (IBLK == 1) THEN
         J = INDEX(NAME,' ')
         IF (INPCI == 0) THEN
            OPEN(UNIT=IK, FILE=NAME(1:J-1)//'.cm', FORM='UNFORMATTED', STATUS=&
               'OLD')
         ELSE
            OPEN(UNIT=IK, FILE=NAME(1:J-1)//'.m', FORM='UNFORMATTED', STATUS=&
               'OLD')
         ENDIF

         READ (IK, IOSTAT=IOS) G92MIX
         IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
            WRITE (6, *)                                              &
                    'File', IK, 'Not a GRASP MIXing Coefficients File;'
            CLOSE(IK)
            STOP
         ENDIF
!
         READ (IK) NELECTOT, NCFTOT, NWTOT, NVECTOT, NVECSIZTOT, NBLOCK1
         IF (NELEC/=NELECTOT .OR. NW/=NWTOT) THEN
!     :    (NCF .NE. NCFT) .OR.
            WRITE (6, *) 'File', IK, 'is not an'
            WRITE (6, *) ' appropriate to Coefficients file'
            CLOSE(IK)
            STOP
         ENDIF
      ENDIF
!
!   Load data from the  rscf.mix  file
!
      WRITE (6, *) 'Loading MIXing Coefficients File ...'
!
      READ (IK) NB, NCF, NVEC, IATJP, IASPA
      CALL ALLOC (EVAL, NVEC, 'EVAL', 'GETMIX')
      CALL ALLOC (EVEC, NCF*NVEC, 'EVEC', 'GETMIX')
      CALL ALLOC (IVEC, NVEC, 'IVEC', 'GETMIX')
      CALL ALLOC (IATJPO, NVEC, 'IATJPO', 'GETMIX')
      CALL ALLOC (IASPAR, NVEC, 'IASPAR', 'GETMIX')
!
!   These arrays are deallocated in mcp
!
      READ (IK) (IVEC(I),I=1,NVEC)
!     READ (IK) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      IATJPO(:NVEC) = IATJP
      IASPAR(:NVEC) = IASPA
      READ (IK) EAV, (EVAL(I),I=1,NVEC)
      READ (IK) ((EVEC(I + (J - 1)*NCF),I=1,NCF),J=1,NVEC)
!
      WRITE (6, *) ' ... load complete;'
!
!   Close the  rscf.mix  file
!
!     CLOSE (IK)
!
      RETURN
      END SUBROUTINE GETMIX
