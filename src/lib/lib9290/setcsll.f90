!***********************************************************************
      SUBROUTINE SETCSLL(NUNIT, NAME, NBLKIN, NBLOCK, NCFBLK, NCFTOT, IDBLK)
!
!  Open, read name file to get nblock, ncfblk(), idblk(), ncftot
!
! Xinghong He 98-06-29
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:34   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!...Modified by Julian Chan 8/12/19
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NUNIT
      INTEGER, INTENT(IN) :: NBLKIN
      INTEGER, INTENT(OUT) :: NBLOCK
      INTEGER, INTENT(OUT) :: NCFTOT
      CHARACTER (LEN = *), INTENT(INOUT) :: NAME
      INTEGER, DIMENSION(*), INTENT(INOUT) :: NCFBLK
      CHARACTER (LEN = 8), DIMENSION(*), INTENT(OUT) :: IDBLK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NCSF, IOS, IERR
      LOGICAL :: FOUND
      CHARACTER :: STR*15, CH*2, LINE3*200, ID*8
!-----------------------------------------------
! Locals

! Look for  <name>
      Print *, 'Entering SETCSLL'

      INQUIRE(FILE=NAME, EXIST=FOUND)
      IF (.NOT.FOUND) THEN
         WRITE (6, *) NAME(1:LEN_TRIM(NAME)), ' does not exist'
         STOP
      ENDIF

! Open it

      CALL OPENFL (NUNIT, NAME, 'FORMATTED', 'OLD', IERR)
      IF (IERR == 1) THEN
         WRITE (6, *) 'Error when opening ', NAME(1:LEN_TRIM(NAME))
         STOP
      ENDIF

! Check the first record of the file; if not as expected, stop

      READ (NUNIT, '(1A15)', IOSTAT=IOS) STR
      IF (IOS/=0 .OR. STR/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         CLOSE(NUNIT)
         STOP
      ENDIF

! Skip next 4 records

      DO I = 1, 4
         READ (NUNIT, *)
      END DO

! Determine the number of blocks in this file

      NBLOCK = 0
      NCSF = 0

      IOS = 0
      DO WHILE(IOS == 0)
         READ (NUNIT, '(A)', IOSTAT=IOS) LINE3
       IF (INDEX(LINE3,'(')/=0) THEN
         READ (NUNIT, '(A)', IOSTAT=IOS) LINE3
       ENDIF
         IF (INDEX(LINE3,'*')/=0 .OR. IOS/=0) THEN
            !.. a new block has been found
            NBLOCK = NBLOCK + 1
            WRITE (6, *) 'Block ', NBLOCK, ',  ncf = ', NCSF
            IF (NBLOCK > NBLKIN) THEN
               WRITE (6, *) 'setcsll: Too many blocks(', NBLOCK, ')'
               WRITE (6, *) 'Maximum allowed is ', NBLKIN
               STOP
            ENDIF
            IDBLK(nblock) = ID
            NCFBLK(NBLOCK) = NCSF
            NCSF = 0
            IF (IOS == 0) CYCLE
         ELSE
       IF(INDEX(LINE3,'+')/=0 .OR. INDEX(LINE3,'-')/=0)THEN
           NCSF = NCSF + 1
                I = LEN_TRIM(LINE3)
                ID = LINE3(I-4:I)
         ENDIF
         ENDIF
      END DO

! Obtain ncftot

      NCFTOT = SUM(NCFBLK(:NBLOCK))

     WRITE(*,*) 'NCFBLK NBLOCK',NCFBLK(1),NCFBLK(2),NCFBLK(3),NBLOCK

      RETURN
      END SUBROUTINE SETCSLL
