!     last edited July 31, 1996
      subroutine lasa1(fil, rad, pop, skal, slut) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reada_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: fil 
      integer  :: skal 
      logical  :: slut 
      character  :: rad*200 
      integer  :: pop(15,0:10,0:1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
!-----------------------------------------------
      if (.not.slut) then 
         do i = 1, 200 
            rad(i:i) = ' ' 
         end do 
         read (fil, 999, end=10) rad 
         call reada (rad, pop, skal, slut) 
         return  
      endif 
   10 continue 
      slut = .TRUE. 
  999 format(a) 
      return  
      end subroutine lasa1 


!     last edited July 30, 1996
      subroutine lasa2(fil, rad2, rad3, stopp, slut) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: fil 
      integer  :: stopp 
      logical , intent(inout) :: slut 
      character  :: rad2*200 
      character  :: rad3*200 
!-----------------------------------------------
      if (.not.slut) then 
         read (fil, 999, end=10) rad2 
!        read(fil,999,end=10) rad2(1:stopp)
         read (fil, 999, end=10) rad3 
!        read(fil,999,end=10) rad3(1:stopp+4)
         return  
      endif 
   10 continue 
      slut = .TRUE. 
  999 format(a) 
      return  
      end subroutine lasa2 


!     last edited July 31, 1996
      subroutine reada(rad1, pop, skal, slut) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(inout) :: skal 
      logical , intent(out) :: slut 
      character , intent(in) :: rad1*200 
      integer , intent(out) :: pop(15,0:10,0:1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, l, antal, stopp 
      character, dimension(0:10) :: orb 
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
      slut = .FALSE. 
      do n = 1, 15 
         pop(n,:min(10,n-1),:1) = 0 
      end do 
      stopp = skal - 1 
      l10: do i = 0, stopp 
         j = 9*i 
         if (rad1(j+3:j+3) == ' ') return  
         skal = i + 1 
         slut = .TRUE. 
         n = ichar(rad1(j+3:j+3)) - ichar('0') 
         if (rad1(j+2:j+2) == '1') n = n + 10 
         if (n<=15 .and. n>=1) then 
            if (rad1(j+7:j+7)==' ' .or. rad1(j+7:j+7)=='0') then 
               do l = 0, min(10,n - 1) 
                  if (rad1(j+4:j+4) /= orb(l)) cycle  
                  slut = .FALSE. 
                  antal = 0 
                  antal = antal + ichar(rad1(j+8:j+8)) - ichar('0') 
                  if (antal > 4*l + 2) then 
                     slut = .TRUE. 
                     return  
                  endif 
                  if (rad1(j+5:j+5)=='-' .or. l==0) then 
                     pop(n,l,0) = antal 
                  else 
                     pop(n,l,1) = antal 
                  endif 
                  cycle  l10 
               end do 
            else 
               do l = 0, min(10,n - 1) 
                  if (rad1(j+4:j+4) /= orb(l)) cycle  
                  slut = .FALSE. 
                  antal = 10*(ichar(rad1(j+7:j+7))-ichar('0')) 
                  antal = antal + ichar(rad1(j+8:j+8)) - ichar('0') 
                  if (antal > 4*l + 2) then 
                     slut = .TRUE. 
                     return  
                  endif 
                  if (rad1(j+5:j+5)=='-' .or. l==0) then 
                     pop(n,l,0) = antal 
                  else 
                     pop(n,l,1) = antal 
                  endif 
                  cycle  l10 
               end do 
            endif 
         else 
            slut = .TRUE. 
            return  
         endif 
      end do l10 
      return  
      end subroutine reada 
