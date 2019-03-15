      subroutine copy7t9
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7
      integer, parameter :: utfil = 9
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character :: rad11*1000
!-----------------------------------------------
      open(fil_1, file='clist.out', status='unknown', position='asis')
      open(unit=utfil, file='fil1.dat', status='unknown', position='asis')
      do while(.TRUE.)
         read (utfil, 999, end=100) rad11
         write (fil_1, 999) trim(rad11)
      end do
  100 continue
      close(utfil, status='delete')
      close(fil_1)
      return
  999 format(a)
      return
      end subroutine copy7t9
