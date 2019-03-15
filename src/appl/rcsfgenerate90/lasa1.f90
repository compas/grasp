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
