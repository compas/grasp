!***********************************************************************
      subroutine starttime(ncount1, progname)
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:49   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncount1
      CHARACTER (LEN = *), INTENT(IN) :: PROGNAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(8) :: nymduhmsm
      integer :: ncount_rate, ncount_max
      character :: chdate*8, chtime*10, chzone*5, msg*80
!-----------------------------------------------
!
! Calls DATE_AND_TIME to get date, time, zone;
!
!
!              !ccyymmdd  hhmmss.sss  Shhmm
!              !Year Month Day Universal Hour Minute Sesond Millisecond
!
!
! For printing

!      write (6, *) '===================================================='
!      write (6, *) '       ', progname, ': Execution Begins ...'
!      write (6, *) '===================================================='

!=======================================================================
!  Get date, time, zone and print
!=======================================================================

!GG      call date_and_time (chdate, chtime, chzone, nymduhmsm)
!GG      write (6, *) 'Date and Time:'
!GG      Print*, '  Date (Yr/Mon/Day): ',                     &
!GG              chdate(1:4),'/',chdate(5:6),'/',chdate(7:8)
!GG      Print*, '  Time (Hr/Min/Sec): ',                     &
!GG              chtime(1:2),'/',chtime(3:4),'/',chtime(5:10)
!GG      Print*, '  Zone: ',chzone

!=======================================================================
!  Start timing - Record the wall clock
!=======================================================================

      call system_clock (ncount1, ncount_rate, ncount_max)

      return
      end subroutine starttime
