!***********************************************************************
!                                                                      *
      SUBROUTINE ANGDATA(NAME,AVAIL,WHICHONE)
!                                                                      *
!   Checks if the angular files name.IOB (one-body) or                 *
!   name.ITB (Two-body) are available                                  *
!                                                                      *
!                                                                      *
!   Modified by C. Naz\'e  Oct. 2011                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = 24), INTENT(IN) :: NAME
      LOGICAL, INTENT(OUT) :: AVAIL
      INTEGER, INTENT(IN)  :: WHICHONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: FOUND
      INTEGER :: J, NF, istat, IWORD
!-----------------------------------------------
!
      J = INDEX(NAME,' ')
      IF (WHICHONE.EQ.1) THEN
          NF = 50
          INQUIRE (FILE = NAME(1:J-1)//'.IOB', EXIST = FOUND)
      ELSEIF (WHICHONE.EQ.2) THEN
          NF = 51
          INQUIRE (FILE = NAME(1:J-1)//'.ITB', EXIST = FOUND)
      ELSE
         print*,'This message should never appear'
      ENDIF

      IF (.NOT.FOUND) THEN
      IF (WHICHONE.EQ.1) PRINT *, ' One-body angular file not available'
      IF (WHICHONE.EQ.2) PRINT *, ' Two-body angular file not available'
        AVAIL = .FALSE.
        RETURN
      ELSE
        IF (WHICHONE.EQ.1) THEN
          OPEN(UNIT=NF,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN'     &
            ,POSITION='APPEND',FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=NF,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN'     &
            ,POSITION='APPEND',FORM='UNFORMATTED')
        ENDIF
        BACKSPACE (UNIT = NF,IOSTAT=istat)
        READ(NF) IWORD
        REWIND(UNIT = NF)
        IF (IWORD.EQ.-1) THEN
          IF (WHICHONE.EQ.1) PRINT *, ' One-body angular file available'
          IF (WHICHONE.EQ.2) PRINT *, ' Two-body angular file available'
          AVAIL = .TRUE.
        ELSE
          AVAIL = .FALSE.
          IF (WHICHONE.EQ.1) PRINT *,' One-body angular file incomplete'
          IF (WHICHONE.EQ.2) PRINT *,' Two-body angular file incomplete'
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ANGDATA
