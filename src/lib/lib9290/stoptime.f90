!***********************************************************************
      subroutine stoptime(ncount1, progname)
!
! Calls DATE_AND_TIME to get date, time, zone;
!
! Things for timing
!              !ccyymmdd  hhmmss.sss  Shhmm
!              !Year Month Day Universal Hour Minute Sesond Millisecond
!
! For printing
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:50   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ncount1
      CHARACTER (LEN = *), INTENT(IN) :: PROGNAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncount2, ncount_rate, ncount_max, nseconds
      integer, dimension(8) :: nymduhmsm
      character :: chdate*8, chtime*10, chzone*5, str2nds*8, msg*80
!-----------------------------------------------

!=======================================================================
!  Get processor info: myid, nprocs, host name; and print
!=======================================================================

!      write (6, *) '===================================================='
!      write (6, *) '       ', progname, ': Execution Finished ...'
!      write (6, *) '===================================================='
      write (6, *)
      write (6, *) 'Wall time:'

      call system_clock (ncount2, ncount_rate, ncount_max)
      ncount2 = ncount2 - ncount1
      nseconds = ncount2/ncount_rate
      write (str2nds, '(I8)') nseconds
      msg = str2nds//' seconds'
      write (6, *) msg(1:len_trim(msg))

      write (6, *)
      write (6, *) 'Finish Date and Time:'

      call date_and_time (chdate, chtime, chzone, nymduhmsm)

      Print*, '  Date (Yr/Mon/Day): ',                     &
              chdate(1:4),'/',chdate(5:6),'/',chdate(7:8)
      Print*, '  Time (Hr/Min/Sec): ',                     &
              chtime(1:2),'/',chtime(3:4),'/',chtime(5:10)
      Print*, '  Zone: ',chzone

      write (6, *)
      write (6, *) progname//': Execution complete.'

      return
      end subroutine stoptime
