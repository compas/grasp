!***********************************************************************
!                                                                      *
      SUBROUTINE CONVRT(INTNUM, CNUM, LENTH)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 22 Sep 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:53   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,   INTENT(IN)    :: INTNUM
      INTEGER,   INTENT(OUT)   :: LENTH
      CHARACTER, INTENT(INOUT) :: CNUM*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=6)                  :: FORM
      CHARACTER(LEN=2), DIMENSION(0:10) :: C1020
      CHARACTER, DIMENSION(9)           :: C19
!
      DATA C19 /'1','2','3','4','5','6','7','8','9'/
      DATA C1020 /'10','11','12','13','14','15','16','17','18','19','20'/
!-----------------------------------------------
      IF (INTNUM < 0) THEN
         LENTH = LOG10(DBLE((-INTNUM))) + 2
      ELSE IF (INTNUM == 0) THEN
         LENTH = 1
      ELSE
         LENTH = LOG10(DBLE(INTNUM)) + 1
      ENDIF
!
!   Ensure that the length of CNUM as dimensioned is adequate;
!   stop with an error message if it isn't
!
      IF (LENTH > LEN(CNUM)) THEN
         WRITE (6, *) 'CONVRT: Length of CNUM inadeuate.'
         ERROR STOP
      ELSE
         IF (LENTH <= 9) THEN
            FORM = '(1I'//C19(LENTH)//')'
            WRITE (CNUM(1:LENTH), FORM(1:5)) INTNUM
         ELSE
            FORM = '(1I'//C1020(LENTH-10)//')'
            WRITE (CNUM(1:LENTH), FORM(1:6)) INTNUM
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE CONVRT
