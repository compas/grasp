!***********************************************************************
!                                                                      *
      SUBROUTINE LDLBL2 (NAME)
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
      CHARACTER, INTENT(IN) :: NAME*24
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
           STATUS='OLD',IOSTAT = IOPEN_STATUS2)
      IF(IOPEN_STATUS2 == 0) THEN
         READ (31,'(1A15)',IOSTAT = IOS) RECORD
         IF (IOS /= 0) THEN
            PRINT *, 'Not a i *.lsj.lbl  File;'
            CLOSE (31)
            STOP
         ELSE
            CALL ALLOC (Lev_POS_2,  NVECTOTF,'Lev_POS_2',  'LDLBL2')
            CALL ALLOC (Lev_J_2,    NVECTOTF,'Lev_J_2',    'LDLBL2')
            CALL ALLOC (Lev_Par_2,  NVECTOTF,'Lev_Par_2',  'LDLBL2')
            CALL ALLOC (RLev_ENER_2,NVECTOTF,'RLev_ENER_2','LDLBL2')
!CPJ        CALL ALLOC (string_CSF2,NVECTOTF,'string_CSF2','LDLBL2')
            allocate(string_CSF2(NVECTOTF))
!
            ICount = 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)     &
              Lev_Pos_2(ICount),Lev_J_2(ICount),Lev_Par_2(ICount),    &
              RLev_ENER_2(ICount)
            IF (IOS .NE. 0) GO TO 1
!
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF2(ICount)
!
    2       READ (31,'(1X,I2)',IOSTAT = IOS),ITEST
            IF (IOS .NE. 0) GO TO 1
            IF (ITEST .EQ. 0) GO TO 2
            BACKSPACE 31
            ICount = ICount + 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)     &
              Lev_Pos_2(ICount),Lev_J_2(ICount),Lev_Par_2(ICount),    &
              RLev_ENER_2(ICount)
            IF (IOS .NE. 0) GO TO 1
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF2(ICount)
            GO TO 2
         ENDIF
      ENDIF
    1 CONTINUE
      CLOSE (31)
      RETURN
      END SUBROUTINE LDLBL2
