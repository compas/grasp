!     last edited July 31, 1996
      subroutine mergeb(antal)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use test_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: antal
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10,0:1) :: pop1, pop2, popo
      integer :: i, j, k
      logical :: p1, p2, slut1, slut2
!-----------------------------------------------
      slut1 = .FALSE.
      slut2 = .FALSE.
      antal = 0
      open(unit=22, status='scratch', position='asis')
      do i = 1, 15
         read (20, 5000, end=2) (pop1(i,j,0),j=0,min(10,i - 1))
         read (20, 5000, end=2) (pop1(i,j,1),j=0,min(10,i - 1))
      end do
      go to 3
    2 continue
      slut1 = .TRUE.
    3 continue
      do i = 1, 15
         read (21, 5000, end=5) (pop2(i,j,0),j=0,min(10,i - 1))
         read (21, 5000, end=5) (pop2(i,j,1),j=0,min(10,i - 1))
      end do
      go to 6
    5 continue
      slut2 = .TRUE.
    6 continue
   10 continue
      if (.not.slut1 .and. .not.slut2) then
         call test (p1, p2, pop1, pop2, 15)
         if (p1) then
            do i = 1, 15
               popo(i,:min(10,i-1),:1) = pop1(i,:min(10,i-1),:1)
            end do
            do i = 1, 15
               write (22, 5000) (pop1(i,j,0),j=0,min(10,i - 1))
               write (22, 5000) (pop1(i,j,1),j=0,min(10,i - 1))
            end do
            do i = 1, 15
               read (20, 5000, end=21) (pop1(i,j,0),j=0,min(10,i - 1))
               read (20, 5000, end=21) (pop1(i,j,1),j=0,min(10,i - 1))
            end do
            go to 22
   21       continue
            slut1 = .TRUE.
   22       continue
            if (p2) then
               do i = 1, 15
                  read (21, 5000, end=23) (pop2(i,j,0),j=0,min(10,i - 1))
                  read (21, 5000, end=23) (pop2(i,j,1),j=0,min(10,i - 1))
               end do
               go to 10
   23          continue
               slut2 = .TRUE.
            endif
         else if (p2) then
            do i = 1, 15
               popo(i,:min(10,i-1),:1) = pop2(i,:min(10,i-1),:1)
            end do
            if (.not.slut2) then
               do i = 1, 15
                  write (22, 5000) (pop2(i,j,0),j=0,min(10,i - 1))
                  write (22, 5000) (pop2(i,j,1),j=0,min(10,i - 1))
               end do
            endif
            do i = 1, 15
               read (21, 5000, end=53) (pop2(i,j,0),j=0,min(10,i - 1))
               read (21, 5000, end=53) (pop2(i,j,1),j=0,min(10,i - 1))
            end do
            go to 10
   53       continue
            slut2 = .TRUE.
         else
            write (*, *) 'fatal error'
            stop
         endif
         go to 10
      else if (.not.slut1 .and. slut2) then
   70    continue
         do i = 1, 15
            write (22, 5000) (pop1(i,j,0),j=0,min(10,i - 1))
            write (22, 5000) (pop1(i,j,1),j=0,min(10,i - 1))
         end do
         do i = 1, 15
            read (20, 5000, end=71) (pop1(i,j,0),j=0,min(10,i - 1))
            read (20, 5000, end=71) (pop1(i,j,1),j=0,min(10,i - 1))
         end do
         go to 70
   71    continue
         slut1 = .TRUE.
      else if (slut1 .and. .not.slut2) then
   80    continue
         do i = 1, 15
            write (22, 5000) (pop2(i,j,0),j=0,min(10,i - 1))
            write (22, 5000) (pop2(i,j,1),j=0,min(10,i - 1))
         end do
         do i = 1, 15
            read (21, 5000, end=81) (pop2(i,j,0),j=0,min(10,i - 1))
            read (21, 5000, end=81) (pop2(i,j,1),j=0,min(10,i - 1))
         end do
         go to 80
   81    continue
         slut2 = .TRUE.
      endif
      rewind (22)
      close(20)
      close(21)
      open(unit=20, status='scratch', position='asis')
  580 continue
      do i = 1, 15
         read (22, 5000, end=999) (pop2(i,j,0),j=0,min(10,i - 1))
         read (22, 5000, end=999) (pop2(i,j,1),j=0,min(10,i - 1))
      end do
      do i = 1, 15
         write (20, 5000) (pop2(i,j,0),j=0,min(10,i - 1))
         write (20, 5000) (pop2(i,j,1),j=0,min(10,i - 1))
      end do
      antal = antal + 1
      go to 580
  999 continue
      close(22)
      rewind (20)
      return
 5000 format(11i2)
      return
      end subroutine mergeb
