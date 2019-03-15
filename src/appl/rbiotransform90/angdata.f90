!***********************************************************************
!                                                                      *
      SUBROUTINE ANGDATA(NAME, AVAIL, KAMAX)
!                                                                      *
!   Checks if the angular file name.T is available and appropriate     *
!                                                                      *
!   Written by Per Jonsson                      6 March 1997           *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE orb_C,           ONLY: NCF, NW
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)      :: KAMAX
      LOGICAL, INTENT(OUT)     :: AVAIL
      CHARACTER, INTENT(INOUT) :: NAME*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, NF, NCFD, NWD, KAMAXD
      LOGICAL :: FOUND
!-----------------------------------------------
!
      J = INDEX(NAME,' ')
      INQUIRE(FILE=NAME(1:J-1)//'.TB', EXIST=FOUND)
      IF (.NOT.FOUND) THEN
         WRITE (6, *) ' Angular file not available'
         AVAIL = .FALSE.
         RETURN
      ELSE
!
!  Open the file and check if it is appropriate for the present case
!
         NF = 200
         OPEN(UNIT=NF, FILE=NAME(1:J-1)//'.TB', STATUS='OLD', FORM=&
            'UNFORMATTED', POSITION='asis')
         REWIND (NF)
         READ (NF) NCFD, NWD, KAMAXD
         IF (.NOT.(NCFD==NCF .AND. NWD==NW .AND. KAMAXD==KAMAX)) THEN
            WRITE (6, *) ' Angular file not appropriate'
            AVAIL = .FALSE.
            RETURN
         ELSE
            WRITE (6, *) ' Angular data read from file'
            AVAIL = .TRUE.
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE ANGDATA
