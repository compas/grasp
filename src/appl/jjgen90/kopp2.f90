!     last edited July 30, 1996
      subroutine kopp2(pos, rad3, j, jprim, par, antko)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: pos
      integer , intent(in) :: par
      character , intent(out) :: rad3*200
      integer , intent(in) :: j(20)
      integer , intent(in) :: jprim(20)
      integer , intent(in) :: antko(20)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, tal
      logical :: first
!-----------------------------------------------

      first = .TRUE.
      do i = 1, 200
         rad3(i:i) = ' '
      end do
      i = max(1,pos - 1)
      k = 9*pos - 2
      if (j(i) == 2*(j(i)/2)) then
         tal = j(i)/20
         if (tal /= 0) rad3(k+2:k+2) = char(48 + tal)
         tal = j(i)/2 - tal*10
         rad3(k+3:k+3) = char(48 + tal)
      else
         tal = j(i)/10
         if (tal /= 0) rad3(k:k) = char(48 + tal)
         tal = j(i) - tal*10
         rad3(k+1:k+1) = char(48 + tal)
         rad3(k+2:k+3) = '/2'
      endif
      if (par == 0) then
         rad3(k+4:k+4) = '+'
      else
         rad3(k+4:k+4) = '-'
      endif
      if (pos > 2) then
         if (jprim(1)/=0 .or. .not.(jprim(1)==0 .and. antko(1)==1)) first = &
            .FALSE.
         do i = 1, pos - 2
            if (first .and. (jprim(i+1)/=0 .or. jprim(i+1)==0 .and. antko(i+1)&
               /=1)) then
               first = .FALSE.
            else if (.not.first .and. jprim(i+1)/=0) then
               k = 9*(i + 1)
               if (j(i) == 2*(j(i)/2)) then
                  tal = j(i)/20
                  if (tal /= 0) rad3(k+2:k+2) = char(48 + tal)
                  tal = j(i)/2 - tal*10
                  rad3(k+3:k+3) = char(48 + tal)
               else
                  tal = j(i)/10
                  if (tal /= 0) rad3(k:k) = char(48 + tal)
                  tal = j(i) - tal*10
                  rad3(k+1:k+1) = char(48 + tal)
                  rad3(k+2:k+3) = '/2'
               endif
            endif
         end do
      endif
      return
      end subroutine kopp2
