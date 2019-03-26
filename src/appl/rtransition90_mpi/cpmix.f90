!************************************************************************
      SUBROUTINE CPMIX(NAME, INPCI)
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: INPCI
      CHARACTER, INTENT(INOUT) :: NAME(2)*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOS, NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK, IBLK, IB, &
         NCF, NVEC, IATJP, IASPA, I, J
      REAL(DOUBLE),DIMENSION(:), pointer :: EVAL, EVEC
      INTEGER,DIMENSION(:), pointer :: IVEC
      REAL(DOUBLE) :: EAV
      CHARACTER :: G92MIX*6
!-----------------------------------------------



      NAME(2) = TRIM(NAME(1))//'_CP'
      IF (INPCI == 0) THEN
         OPEN(UNIT=78, FILE=TRIM(NAME(2))//'.cbm', FORM='UNFORMATTED', STATUS=&
            'UNKNOWN', POSITION='asis')
      ELSE
         OPEN(UNIT=78, FILE=TRIM(NAME(2))//'.bm', FORM='UNFORMATTED', STATUS=&
            'UNKNOWN', POSITION='asis')
      ENDIF
      IF (INPCI == 0) THEN
         OPEN(UNIT=68, FILE=TRIM(NAME(1))//'.cbm', FORM='UNFORMATTED', STATUS=&
            'OLD', POSITION='asis')
      ELSE
         OPEN(UNIT=68, FILE=TRIM(NAME(1))//'.bm', FORM='UNFORMATTED', STATUS=&
            'OLD', POSITION='asis')
      ENDIF
      READ (68, IOSTAT=IOS) G92MIX
      IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
         WRITE (*, *) 'Not a GRASP mixing file'
         STOP
      ENDIF
      WRITE (78) G92MIX
      READ (68) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK
      WRITE (78) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK

      DO IBLK = 1, NBLOCK
         READ (68) IB, NCF, NVEC, IATJP, IASPA
         WRITE (78) IB, NCF, NVEC, IATJP, IASPA
         CALL ALLOC (EVAL, NVEC, 'EVAL', 'CPMIX')
         CALL ALLOC (EVEC, NCF*NVEC, 'EVEC', 'CPMIX')
         CALL ALLOC (IVEC, NVEC, 'IVEC', 'CPMIX')
         READ (68) (IVEC(I),I=1,NVEC)
         READ (68) EAV, (EVAL(I),I=1,NVEC)
         READ (68) ((EVEC(I + (J - 1)*NCF),I=1,NCF),J=1,NVEC)

         WRITE (78) (IVEC(I),I=1,NVEC)
         WRITE (78) EAV, (EVAL(I),I=1,NVEC)
         WRITE (78) ((EVEC(I + (J - 1)*NCF),I=1,NCF),J=1,NVEC)
         CALL DALLOC (EVAL, 'EVAL', 'CPMIX')
         CALL DALLOC (EVEC, 'EVEC', 'CPMIX')
         CALL DALLOC (IVEC, 'IVEC', 'CPMIX')
      END DO
      NAME(2) = NAME(1)
      CLOSE(68)
      CLOSE(78)
      RETURN
      END SUBROUTINE CPMIX
