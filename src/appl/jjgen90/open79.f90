      subroutine open79(i)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7
      integer, parameter :: utfil = 9
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i1
!-----------------------------------------------
      i1 = mod(i,2)
      if (i1 == 0) then
         open(fil_1, file='fil1.dat', status='unknown', position='asis')
         open(unit=utfil, file='clist.out', status='unknown', position='asis')
      else
         open(fil_1, file='clist.out', status='unknown', position='asis')
         open(unit=utfil, file='fil1.dat', status='unknown', position='asis')
      endif
      close(utfil)
      rewind (fil_1)
      return
      end subroutine open79
