!***********************************************************************
!                                                                      *
      SUBROUTINE SETRES(ISOFILE, RWFFILE, IDBLK) 
!                                                                      *
!   Open, check, load data from the  .res  file.                       *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified by Xinghong                  Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE memory_man
      USE default_C
      USE where_C
      USE hblock_C
      USE iccu_C
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE openfl_I 
      USE lodres_I 
      USE getcid_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: ISOFILE*(*) 
      CHARACTER  :: RWFFILE*(*) 
      CHARACTER(LEN=8), DIMENSION(*)  :: IDBLK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*11, PARAMETER :: FORM = 'UNFORMATTED' 
      CHARACTER*7, PARAMETER :: STATUS = 'UNKNOWN' 
      CHARACTER*6, PARAMETER :: RESTITLE = 'R92RES' 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IERR, IOS 
      LOGICAL :: FOUND, RESTRT 
      CHARACTER :: R92RES*6, DEFNAM*11, IDSTRING*3 
!-----------------------------------------------
!
! Compose the "rci.res" file name
!
      DEFNAM = 'rci.res' 
!
! Ask if this is a restart
!
      IF (NDEF /= 0) THEN 
         WRITE (ISTDE, *) 'Restarting RCI90 ?' 
         RESTRT = GETYN() 
      ELSE 
         RESTRT = .FALSE. 
      ENDIF 
!      IF (RESTRT) THEN
!        WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
!      ELSE
!        WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
!      END IF
!
! Do some settings and checks
!
      IF (RESTRT) THEN 
!         ...Restart, make sure file exist
         INQUIRE(FILE=DEFNAM, EXIST=FOUND) 
         IF (.NOT.FOUND) STOP 'setres: .res does not exist' 
      ENDIF 
!
! Open the .res file
!
      CALL OPENFL (IMCDF, DEFNAM, FORM, STATUS, IERR) 
      IF (IERR /= 0) STOP 'setres: Error openning .res file' 
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      CALL ALLOC (ICCUTBLK, NBLOCK, 'ICCUTBLK', 'SETRES') 
 
      IF (RESTRT) THEN 
!        ...Check the signature of the file
         READ (IMCDF, IOSTAT=IOS) R92RES 
         IF (IOS/=0 .OR. R92RES/=RESTITLE) THEN 
            CLOSE(IMCDF) 
            STOP 'setres: Not RCI92 .res file' 
         ENDIF 
 
!         ...Read and check restart information
         CALL LODRES 
 
      ELSE 
 
!         ...Write the file header
!         ...Generate the first part of the .res file
         WRITE (IMCDF) RESTITLE 
         CALL GETCID (ISOFILE, RWFFILE, IDBLK) 
 
      ENDIF 
 
      RETURN  
      END SUBROUTINE SETRES 
