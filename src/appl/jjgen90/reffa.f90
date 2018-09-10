!     last edited July 31, 1996
      subroutine reffa(posn, posl) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(inout) :: posn(110) 
      integer , intent(inout) :: posl(110) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: logfil = 31 
      integer, parameter :: reffil = 18 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10) :: stat 
      integer :: i, n, l, num 
      logical :: ok 
      character , dimension(0:10) :: orb 
      character :: m, x, y 
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
 
      do n = 1, 15 
         stat(n,:min(n-1,10)) = 0 
      end do 
      write (*, 200) 'Default, reverse, symmetry or user specified ordering?', &
         ' (*/r/s/u)' 
      read (*, 1000) x 
      if (x=='u' .or. x=='U') then 
         write (logfil, *) 'User specified ordering.' 
         inquire(file='clist.ref', exist=ok) 
         if (ok) open(unit=reffil, status='old', file='clist.ref', position=&
            'asis') 
         l = -1 
         num = 1 
         if (.not.ok) then 
            write (*, 200) 'No reference file present! ', &
               'The couplings will appear in standard order.' 
         else 
            write (*, 200) 'Reference file present!' 
   20       continue 
            read (reffil, 1000, end=40) m, x, y 
            n = ichar(m) - ichar('0') 
            if (x>='0' .and. x<='9') then 
               n = n*10 + ichar(x) - ichar('0') 
               x = y 
            endif 
            do i = 0, 10 
               if (orb(i) /= x) cycle  
               l = i 
            end do 
            if (l==(-1) .or. n<0 .or. n>15 .or. n<=l .or. l>10) go to 40 
            if (stat(n,l) /= 0) then 
               write (*, 200) 'The same orbital appeared more than once!' 
               l = -1 
               go to 20 
            endif 
            posn(num) = n 
            posl(num) = l 
            stat(n,l) = num 
            num = num + 1 
            l = -1 
            go to 20 
   40       continue 
            if (num == 1) then 
               write (*, 200) 'The program failed reading the order of ', &
                  'the customized coupling scheme.' 
            else 
               write (*, 200) 'The couplings will ', &
                  'be made in the following customized order:' 
               if (num == 2) then 
                  write (*, 100) posn(1), orb(posl(1)) 
               else 
                  write (*, 100) posn(1), orb(posl(1)), (',',posn(i),orb(posl(i&
                     )),i=2,num - 1) 
               endif 
            endif 
         endif 
         do n = 1, 15 
            do l = 0, min(n - 1,10) 
               if (stat(n,l) /= 0) cycle  
               posn(num) = n 
               posl(num) = l 
               num = num + 1 
            end do 
         end do 
         close(reffil) 
         write (*, 200) 
      else if (x=='s' .or. x=='S') then 
         write (logfil, *) 'Symmetry ordering.' 
         num = 1 
         do l = 0, 10 
            do n = l + 1, 15 
               posn(num) = n 
               posl(num) = l 
               num = num + 1 
            end do 
         end do 
      else if (x=='r' .or. x=='R') then 
         write (logfil, *) 'Reverse ordering.' 
         num = 1 
         do n = 15, 1, -1 
            do l = min(n - 1,10), 0, -1 
               posn(num) = n 
               posl(num) = l 
               num = num + 1 
            end do 
         end do 
      else 
         write (logfil, *) 'Standard ordering.' 
         num = 1 
         do n = 1, 15 
            do l = 0, min(n - 1,10) 
               posn(num) = n 
               posl(num) = l 
               num = num + 1 
            end do 
         end do 
      endif 
      return  
  100 format(' ',110(i2,a,a)) 
  200 format(' ',2a) 
 1000 format(3a) 
      return  
      end subroutine reffa 
