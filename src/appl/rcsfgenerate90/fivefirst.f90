!     last edited Februar 20, 1996
      subroutine fivefirst(slut1, slut2, posn, posl)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical , intent(inout) :: slut1
      logical , intent(inout) :: slut2
      integer , intent(in) :: posn(110)
      integer , intent(in) :: posl(110)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j, n, l, pos, stopp
      logical, dimension(15,0:10,0:1) :: med
      character :: rad0*1000, rad1*1000, rad2*1000
      character, dimension(0:10) :: orb
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/
!-----------------------------------------------


      read (7, 999, end=100)
      write (9, 999) 'Core subshells:'
      read (7, 999, end=101) rad0
      stopp = 0
      do i = 0, 210
         if (ichar(rad0(stopp+3:stopp+3))>=ichar('0') .and. ichar(rad0(stopp+3:&
            stopp+3))<=ichar('9')) then
            stopp = stopp + 5
         else
            exit
         endif
      end do
      if (stopp /= 0) then
         write (9, 999) rad0(1:stopp)
      else
         write (9, 999)
      endif
      write (9, 999) 'Peel subshells:'
      read (7, 999, end=102)
      read (7, 999, end=102) rad1
      read (7, 999, end=102)
      if (.not.slut2) then
         read (8, 999, end=90)
         read (8, 999, end=90)
         read (8, 999, end=90)
         read (8, 999, end=90) rad2
         read (8, 999, end=90)
         do i = 1, 15
            med(i,:min(10,i-1),0) = .FALSE.
            med(i,:min(10,i-1),1) = .FALSE.
         end do
         do i = 1, 205
            pos = 5*i
            n = ichar(rad1(pos-2:pos-2)) - ichar('0')
            if (rad1(pos-3:pos-3) == '1') n = n + 10
            l = -1
            if (n>=1 .and. n<=15) then
               do j = 0, min(10,n - 1)
                  if (rad1(pos-1:pos-1) /= orb(j)) cycle
                  l = j
               end do
            endif
            if (l == (-1)) exit
            if (rad1(pos:pos)=='-' .or. l==0) then
               med(n,l,0) = .TRUE.
            else
               med(n,l,1) = .TRUE.
            endif
         end do
         do i = 1, 205
            pos = 5*i
            n = ichar(rad2(pos-2:pos-2)) - ichar('0')
            if (rad2(pos-3:pos-3) == '1') n = n + 10
            l = -1
            if (n>=1 .and. n<=15) then
               do j = 0, min(10,n - 1)
                  if (rad2(pos-1:pos-1) /= orb(j)) cycle
                  l = j
               end do
            endif
            if (l == (-1)) exit
            if (rad2(pos:pos)=='-' .or. l==0) then
               med(n,l,0) = .TRUE.
            else
               med(n,l,1) = .TRUE.
            endif
         end do
         pos = 3
         do k = 1, 110
            i = posn(k)
            j = posl(k)
            if (med(i,j,0)) then
               rad0(pos-2:pos+2) = '     '
               if (i < 10) then
                  rad0(pos:pos) = char(i + ichar('0'))
               else
                  rad0(pos:pos) = char(i + ichar('0') - 10)
                  rad0(pos-1:pos-1) = '1'
               endif
               rad0(pos+1:pos+1) = orb(j)
               if (j /= 0) rad0(pos+2:pos+2) = '-'
               pos = pos + 5
            endif
            if (.not.med(i,j,1)) cycle
            rad0(pos-2:pos+2) = '     '
            if (i < 10) then
               rad0(pos:pos) = char(i + ichar('0'))
            else
               rad0(pos:pos) = char(i + ichar('0') - 10)
               rad0(pos-1:pos-1) = '1'
            endif
            rad0(pos+1:pos+1) = orb(j)
            pos = pos + 5
         end do
         write (9, 999) rad0(1:pos-3)
         write (9, 999) 'CSF(s):'
         return
      endif
   90 continue
      slut2 = .TRUE.
      stopp = 0
      do i = 0, 210
         if (ichar(rad1(stopp+3:stopp+3))>=ichar('0') .and. ichar(rad1(stopp+3:&
            stopp+3))<=ichar('9')) then
            stopp = stopp + 5
         else
            exit
         endif
      end do
      if (stopp /= 0) then
         write (9, 999) rad1(1:stopp)
      else
         write (9, 999)
      endif
      write (9, 999) 'CSF(s):'
      return
  100 continue
      slut1 = .TRUE.
      write (9, 999) 'Core subshells:'
      read (8, 999, end=200) rad0
  101 continue
      if (.not.slut1) then
         slut1 = .TRUE.
         read (8, 999, end=200) rad0
      endif
      stopp = 0
      do i = 0, 210
         if (ichar(rad0(stopp+3:stopp+3))>=ichar('0') .and. ichar(rad0(stopp+3:&
            stopp+3))<=ichar('9')) then
            stopp = stopp + 5
         else
            exit
         endif
      end do
      if (stopp /= 0) then
!PJ         write(9,999) rad1(1:stopp)
         do i = 1,1000
           if (rad1(i:i).eq.':') rad1(i-1:i) = '10'
           if (rad1(i:i).eq.';') rad1(i-1:i) = '11'
           if (rad1(i:i).eq.'<') rad1(i-1:i) = '12'
           if (rad1(i:i).eq.'=') rad1(i-1:i) = '13'
           if (rad1(i:i).eq.'>') rad1(i-1:i) = '14'
           if (rad1(i:i).eq.'?') rad1(i-1:i) = '15'
         end do
         write(9,999) trim(rad1)
!PJ

!        write (9, 999) rad0(1:stopp)
      else
         write (9, 999)
      endif
      write (9, 999) 'Peel subshells:'
  102 continue
      if (.not.slut1) then
         slut1 = .TRUE.
         read (8, 999, end=200)
         read (8, 999, end=200) rad2
      endif
      stopp = 0
      do i = 0, 210
         if (ichar(rad2(stopp+3:stopp+3))>=ichar('0') .and. ichar(rad2(stopp+3:&
            stopp+3))<=ichar('9')) then
            stopp = stopp + 5
         else
            exit
         endif
      end do
      if (stopp /= 0) then
         write (9, 999) rad2(1:stopp)
      else
         write (9, 999)
      endif
      write (9, 999) 'CSF(s):'
      return
  200 continue
      slut2 = .TRUE.
      return
  999 format(a)
      return
      end subroutine fivefirst
