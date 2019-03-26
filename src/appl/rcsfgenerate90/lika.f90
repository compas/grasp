
!     last edited September 23, 1995
      logical function lika (pop0, pop1)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: pop0(15,0:10,0:1)
      integer , intent(in) :: pop1(15,0:10,0:1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k
      logical :: dum
!-----------------------------------------------
      dum = .TRUE.
      l10: do i = 1, 15
         do j = 0, min(10,i - 1)
            do k = 0, 1
               dum = dum .and. pop0(i,j,k)==pop1(i,j,k)
               if (dum) cycle
               exit  l10
            end do
         end do
      end do l10
      lika = dum
      return
      end function lika
