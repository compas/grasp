!     last edited July 30, 1996
      subroutine kopp1(pos, rad2, j, s, antko)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: pos
      character , intent(out) :: rad2*200
      integer , intent(in) :: j(20)
      integer , intent(in) :: s(20)
      integer , intent(in) :: antko(20)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, tal
!-----------------------------------------------
      do i = 1, 200
         rad2(i:i) = ' '
      end do
      do i = 1, pos
         k = 9*i
         if (j(i) == 2*(j(i)/2)) then
            if (.not.(j(i)==0 .and. antko(i)==1)) then
               if (s(i) /= (-1)) then
                  rad2(k-4:k-4) = ';'
                  rad2(k-5:k-5) = char(48 + s(i))
               endif
               tal = j(i)/20
               if (tal /= 0) then
!GG                  rad2(k:k) = char(48 + tal)
                  rad2(k-1:k-1) = char(48+tal)
               end if
               tal = j(i)/2 - tal*10
               rad2(k:k) = char(48 + tal)
            endif
         else
            if (s(i) /= (-1)) then
               rad2(k-4:k-4) = ';'
               rad2(k-5:k-5) = char(48 + s(i))
            endif
            tal = j(i)/10
            if (tal /= 0) rad2(k-3:k-3) = char(48 + tal)
            tal = j(i) - tal*10
            rad2(k-2:k-2) = char(48 + tal)
            rad2(k-1:k) = '/2'
         endif
      end do
      return
      end subroutine kopp1
