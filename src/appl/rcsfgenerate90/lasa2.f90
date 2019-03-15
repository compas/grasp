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
