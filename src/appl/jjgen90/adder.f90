!     last edited July 31, 1996
      subroutine adder(closed, med, slut, anel, par, expand)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lockad_I
      use lasa1_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: anel
      integer , intent(out) :: par
      logical  :: slut
      logical  :: expand
      logical  :: closed(15,0:10)
      logical  :: med(15,0:10)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7
      integer, parameter :: fil_2 = 8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10,0:1) :: pop
      integer :: skal, i, j, kl, nr
      logical :: finns
      character :: rad1*500, rad2*500, rad3*500
!-----------------------------------------------
      skal = 20
      inquire(file='clist.inp', exist=finns)
      if (finns) then
         if (.not.expand) then
            open(unit=fil_1, file='clist.inp', status='old', position='asis')
         else
            open(unit=fil_2, file='clist.inp', status='old', position='asis')
         endif
         slut = .FALSE.
         call lockad (closed, med, slut, expand)
         if (.not.slut) then
            if (expand) then
               read (fil_2, *, end=99)
               call lasa1 (fil_2, rad1, pop, skal, slut)
            else
               read (fil_1, *, end=99)
               call lasa1 (fil_1, rad1, pop, skal, slut)
            endif
         endif
         if (.not.slut) then
            anel = 0
            par = 0
            do i = 1, 15
               do j = 0, min(10,i - 1)
                  if (closed(i,j)) then
                     anel = anel + 2 + 4*j
                  else
                     anel = anel + pop(i,j,0) + pop(i,j,1)
                     par = mod(par + j*(pop(i,j,0)+pop(i,j,1)),2)
                  endif
               end do
            end do
            if (expand) then
               read (fil_2, 100, end=99) rad2
               read (fil_2, 100, end=99) rad3
            else
               read (fil_1, 100, end=99) rad2
               read (fil_1, 100, end=99) rad3
            endif
            kl = skal*9
            if (rad3(kl:kl) /= '/') then
               if (rad3(kl:kl) /= ' ') then
                  nr = 10*(ichar(rad3(kl:kl))-ichar('0'))
               else
                  nr = 0
               endif
               kl = kl + 1
               if (rad3(kl:kl) /= ' ') nr = nr + (ichar(rad3(kl:kl))-ichar('0')&
                  )
            else
               kl = skal*9 - 2
               if (rad3(kl:kl) /= ' ') then
                  nr = 10*(ichar(rad3(kl:kl))-ichar('0'))
               else
                  nr = 0
               endif
               kl = kl + 1
               nr = nr + ichar(rad3(kl:kl)) - ichar('0')
            endif
         endif
         if (expand) then
            rewind (fil_2)
         else
            rewind (fil_1)
         endif
      else
         slut = .TRUE.
      endif
      return
   99 continue
      slut = .TRUE.
      if (expand) then
         close(fil_2)
      else
         close(fil_1)
      endif
      return
  100 format(a)
      return
      end subroutine adder
