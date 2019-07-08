!***********************************************************************
!                                                                      *
      SUBROUTINE LDLBL1 (NAME)
!                                                                      *
!   Open, check and load data from the  .lsj.lbl   file of the         *
!   inital state.                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC.                                        *
!                                                                      *
!   Written by G. Gaigalas,                                            *
!   NIST                                                  May 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man
      USE jj2lsjbio_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER, INTENT(IN) :: NAME*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: WEIGHTS
      CHARACTER    :: RECORD*15
      INTEGER      :: ICount, IOS, ITEST, J
!-----------------------------------------------
!
!   Check the first record of the file; if not as expected, try again
!
      J = INDEX(NAME,' ')
      OPEN (UNIT = 31,FILE=NAME(1:J-1)//'.lsj.lbl',FORM='FORMATTED',  &
           STATUS='OLD',IOSTAT = IOPEN_STATUS1)
      IF(IOPEN_STATUS1 == 0) THEN
         READ (31,'(1A15)',IOSTAT = IOS) RECORD
         IF (IOS /= 0) THEN
            PRINT *, 'Not a i *.lsj.lbl  File;'
            CLOSE (31)
            STOP
         ELSE
            CALL ALLOC (Lev_POS_1,  NVECTOTI,'Lev_POS_1',  'LDLBL1')
            CALL ALLOC (Lev_J_1,    NVECTOTI,'Lev_J_1',    'LDLBL1')
            CALL ALLOC (Lev_Par_1,  NVECTOTI,'Lev_Par_1',  'LDLBL1')
            CALL ALLOC (RLev_ENER_1,NVECTOTI,'RLev_ENER_1','LDLBL1')
!CPJ        CALL ALLOC (string_CSF1,NVECTOTI,'string_CSF1','LDLBL1')
            allocate(string_CSF1(1:NVECTOTI))
!

            ICount = 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)     &
              Lev_Pos_1(ICount),Lev_J_1(ICount),Lev_Par_1(ICount),    &
              RLev_ENER_1(ICount)
            IF (IOS .NE. 0) GO TO 1
!
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF1(ICount)
!
    2       READ (31,'(1X,I2)',IOSTAT = IOS) ITEST
            IF (IOS .NE. 0) GO TO 1
            IF (ITEST .EQ. 0) GO TO 2
            BACKSPACE 31
            ICount = ICount + 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)     &
              Lev_Pos_1(ICount),Lev_J_1(ICount),Lev_Par_1(ICount),    &
              RLev_ENER_1(ICount)
            IF (IOS .NE. 0) GO TO 1
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF1(ICount)
            GO TO 2
         ENDIF
      ENDIF
    1 CONTINUE
      CLOSE (31)
      RETURN
      END SUBROUTINE LDLBL1
