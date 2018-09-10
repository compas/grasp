!***********************************************************************
!                                                                      *
      SUBROUTINE ANGDATA(NAME, AVAIL, JKP, NFILE2) 
!                                                                      *
!   Checks if the angular file name(1).name(2).T is available          *
!   and appropriate                                                    *
!                                                                      *
!   Written by Per Jonsson                      6 March 1997           *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE orb_C, ONLY: NCF, NW, IQA
      USE osc_C, ONLY: NKP, KP
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JKP 
      INTEGER, INTENT(IN) :: NFILE2 
      LOGICAL, INTENT(OUT) :: AVAIL 
      CHARACTER, INTENT(INOUT) :: NAME(2)*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1, J2, IBLKI, IBLKF, NWD, NKPD 
      LOGICAL :: FOUND 
      CHARACTER, DIMENSION(-9:9) :: S*2 
!-----------------------------------------------
 
      S((-9)) = '-9' 
      S((-8)) = '-8' 
      S((-7)) = '-7' 
      S((-6)) = '-6' 
      S((-5)) = '-5' 
      S((-4)) = '-4' 
      S((-3)) = '-3' 
      S((-2)) = '-2' 
      S((-1)) = '-1' 
      S(0) = '+0' 
      S(1) = '+1' 
      S(2) = '+2' 
      S(3) = '+3' 
      S(4) = '+4' 
      S(5) = '+5' 
      S(6) = '+6' 
      S(7) = '+7' 
      S(8) = '+8' 
      S(9) = '+9' 
 
      J1 = INDEX(NAME(1),' ') 
      J2 = INDEX(NAME(2),' ') 
      INQUIRE(FILE=NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'//S(KP(JKP))//'T'&
         , EXIST=FOUND) 
      IF (.NOT.FOUND) THEN 
         WRITE (6, *) 
         WRITE (6, *) ' Angular file not available' 
         AVAIL = .FALSE. 
         RETURN  
      ELSE 
!
!  Open the file and check if it is appropriate for the present case
!
         OPEN(UNIT=NFILE2, FILE=NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'//S(&
            KP(JKP))//'T', STATUS='OLD', FORM='UNFORMATTED', POSITION='asis') 
         REWIND (NFILE2) 
         READ (NFILE2) IBLKI, IBLKF, NWD, NKPD 
         IF (.NOT.(NWD==NW .AND. NKPD==NKP)) THEN 
            WRITE (6, *) ' Angular file not appropriate' 
            AVAIL = .FALSE. 
            CLOSE(NFILE2, STATUS='DELETE') 
            RETURN  
         ELSE 
            REWIND (NFILE2) 
            WRITE (6, *) ' Angular data read from file' 
            AVAIL = .TRUE. 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE ANGDATA 
