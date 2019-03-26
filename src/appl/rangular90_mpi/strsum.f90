!***********************************************************************
!                                                                      *
      SUBROUTINE STRSUM
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      Use hblock_C
      USE def_C
      USE iccu_C
      USE mcp_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE calen_I
      USE convrt_I
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENTH, I
      Character :: CTIME*10, CDATE*8
      CHARACTER :: RECORD*256, CDATA*26
!-----------------------------------------------
!
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
      CALL CALEN (CTIME, CDATE)
      WRITE (24, *) 'GENMCP run at ', CTIME, ' on ', CDATE, '.'
!
!   Write out the basic dimensions of the electron cloud
!
      WRITE (24, *)
      CALL CONVRT (NELEC, RECORD, LENTH)
      WRITE (24, *) 'There are '//RECORD(1:LENTH)//' electrons in the cloud'
      CALL CONVRT (NCF, RECORD, LENTH)
      WRITE (24, *) ' in '//RECORD(1:LENTH)//' relativistic CSFs'
      CALL CONVRT (NW, RECORD, LENTH)
      WRITE (24, *) ' based on '//RECORD(1:LENTH)//' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      IF (DIAG) THEN
         WRITE (24, *)
         WRITE (24, *) 'Only diagonal matrix elements are computed.'
      ELSE
         IF (LFORDR) THEN
            DO I = 1, NBLOCK
               WRITE (24, *)
               CALL CONVRT (ICCUT(I), RECORD, LENTH)
               WRITE (24, *) 'CSFs 1--'//RECORD(1:LENTH)//' constitute'//&
                  ' the zero-order space.'
            END DO
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE STRSUM
