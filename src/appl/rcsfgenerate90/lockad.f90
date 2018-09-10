!     last edited Februar 20, 1996
      subroutine lockad(closed, med, slut, expand) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical , intent(out) :: slut 
      logical , intent(in) :: expand 
      logical , intent(out) :: closed(15,0:10) 
      logical , intent(out) :: med(15,0:10) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7 
      integer, parameter :: fil_2 = 8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, n, l 
      character :: rad*1000 
      character, dimension(0:10) :: orb 
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
      if (expand) then 
         read (fil_2, *, end=40) 
         read (fil_2, 100, end=40) rad 
      else 
         read (fil_1, *, end=40) 
         read (fil_1, 100, end=40) rad 
      endif 
      do n = 1, 15 
         closed(n,1:min(10,n-1)) = .FALSE. 
      end do 
      l30: do i = 0, 205 
         j = i*5 
         n = ichar(rad(j+3:j+3)) - ichar('0') 
         if (n>=1 .and. n<=15) then 
            do l = 0, min(10,n - 1) 
               if (rad(j+4:j+4) /= orb(l)) cycle  
               closed(n,l) = .TRUE. 
               cycle  l30 
            end do 
         else 
            exit  l30 
         endif 
      end do l30 
      if (expand) then 
         read (fil_2, *, end=40) 
         read (fil_2, 100, end=40) rad 
      else 
         read (fil_1, *, end=40) 
         read (fil_1, 100, end=40) rad 
      endif 
      do n = 1, 15 
         med(n,1:min(10,n-1)) = .FALSE. 
      end do 
      l130: do i = 0, 205 
         j = i*5 
         n = ichar(rad(j+3:j+3)) - ichar('0') 
         if (n>=1 .and. n<=15) then 
            do l = 0, min(10,n - 1) 
               if (rad(j+4:j+4) /= orb(l)) cycle  
               med(n,l) = .TRUE. 
               cycle  l130 
            end do 
         else 
            return  
         endif 
      end do l130 
 
      return  
   40 continue 
      slut = .TRUE. 
      return  
  100 format(a) 
      return  
      end subroutine lockad 
