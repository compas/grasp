!     last edited July 30, 1996
      subroutine fivelines(org, locked, closed, first, posn, posl)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical , intent(in) :: first
      integer , intent(inout) :: org(15,0:10)
      integer , intent(in) :: posn(110)
      integer , intent(in) :: posl(110)
      logical , intent(in) :: locked(15,0:10)
      logical , intent(in) :: closed(15,0:10)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j, start, stopp
      character :: rad*1000
      character, dimension(0:10) :: orb
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/
      if (.not.first) then
!        open(unit=8,file='fil2.dat',status='unknown')
         open(unit=8, status='scratch', position='asis')
         write (8, 999)
         write (8, 999)
         write (8, 999)
      else
         open(unit=7, file='fil1.dat', status='unknown', position='asis')
!        open(unit=7,status='scratch')
         write (7, 999) 'Core subshells:'
         do i = 1, 1000
            rad(i:i) = ' '
         end do
         start = -2
         stopp = 0
         do k = 1, 110
            i = posn(k)
            j = posl(k)
            if (.not.closed(i,j)) cycle
            start = start + 5
            stopp = stopp + 5
            if (i < 10) then
               rad(start:start) = char(ichar('0') + i)
            else
               rad(start-1:start-1) = '1'
               rad(start:start) = char(ichar('0') + i - 10)
            endif
            rad(start+1:start+1) = orb(j)
            org(i,j) = 0
            if (j < 1) cycle
            rad(start+2:start+2) = '-'
            start = start + 5
            stopp = stopp + 5
            if (i < 10) then
               rad(start:start) = char(ichar('0') + i)
            else
               rad(start-1:start-1) = '1'
               rad(start:start) = char(ichar('0') + i - 10)
            endif
            rad(start+1:start+1) = orb(j)
         end do
         if (stopp == 0) then
            write (7, 999)
         else
            write (7, 999) rad(1:stopp)
         endif
         write (7, 999) 'Peel subshells:'
      endif
      do i = 1, 1000
         rad(i:i) = ' '
      end do
      start = -2
      stopp = 0
      do k = 1, 110
         i = posn(k)
         j = posl(k)
         if (.not.(.not.(org(i,j)==0 .and. locked(i,j)) .and. .not.closed(i,j))&
            ) cycle
         start = start + 5
         stopp = stopp + 5
         if (i < 10) then
            rad(start:start) = char(ichar('0') + i)
         else
            rad(start-1:start-1) = '1'
            rad(start:start) = char(ichar('0') + i - 10)
         endif
         rad(start+1:start+1) = orb(j)
         if (j < 1) cycle
         rad(start+2:start+2) = '-'
         start = start + 5
         stopp = stopp + 5
         if (i < 10) then
            rad(start:start) = char(ichar('0') + i)
         else
            rad(start-1:start-1) = '1'
            rad(start:start) = char(ichar('0') + i - 10)
         endif
         rad(start+1:start+1) = orb(j)
!           write(*,*) i,rad(1:100)
      end do
      if (first) then
         if (stopp == 0) then
            write (7, 999)
         else
            write (7, 999) rad(1:stopp)
         endif
         write (7, 999) 'CSF(s):'
      else
         if (stopp == 0) then
            write (8, 999)
         else
            write (8, 999) rad(1:stopp)
         endif
         write (8, 999) 'CSF(s):'
      endif
  999 format(a)
      return
      end subroutine fivelines
