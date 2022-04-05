PROGRAM rwfntotxt

   IMPLICIT NONE

   integer, parameter                   :: npts0 = 10000 ! max number of gridpoints
   integer, parameter                   :: nnnw  = 127   ! max no of orbitals
   character(len=1), dimension(10), parameter :: l_str = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm']

   real(kind=8), dimension(npts0, nnnw) :: pg, qg, rg, pgg
   real(kind=8)                         :: energy, a0
   character(len=5), dimension(nnnw)    :: nl, selected_orbs
   character                            :: title*6, orbl*4, n_str*2, new*3, xa*3, leg*50, format_string*50
   integer(kind=4)                      :: np, lp, jp, nn, laky, ll, jj, npts, j, i, nwf, i_rg_max

   write (*, *)
   write (*, *) ' RWFNTOTXT - Jon Grumer, Uppsala University, 2022'
   write (*, *) ' Converts binary GRASP wavefunctions to ASCII file format'
   write (*, *)
   write (*, *) '    Input: rwfn.inp      Output: rwfn.plot'
   write (*, *)
   write (*, *) '    Note: rwfnpyplot can be used to conveniently plot orbitals from the .plot file.'
   write (*, *) '          E.g. >> rwfnpyplot rwfn.plot 2p- 50 y'
   write (*, *) '          plots the large P(2p-) component to screen up to r = 50 au.'
   write (*, *)

   ! optional:
   ! selected_orbs(1:2) = ['6s', '7s']

   open (3, file='rwfn.inp', status='old', form='unformatted')
   open (4, file='rwfn.plot', status='unknown')

   ! write terminal info header
   write(*, '(a32, a9)') 'max ', 'max'
   write(*, '(a5, a7, a14, a7, a10, a6)') 'i', 'nl  ', 'E[a.u]      ', 'grid-n', 'r[a.u]', 'a0'

   ! read binary orbital file
   read (3) title
   if (title .ne. 'G92RWF') then
      print *, 'title = ', title, 'does not match G92RWF'
      stop
   end if

   i = 1
   do
      read (3, end=20) nn, laky, energy, npts

      write (n_str, '(i2)') nn
      if (laky .gt. 0) then
         ll = laky
         jj = -1
      elseif (laky .le. -1) then
         ll = -laky - 1
         jj = 1
      else
         write (*, *) 'Unexpected case when reading GRASP orbital wavefunction file'
         stop
      end if

      ! define nl+/- in orbital array
      if (jj .eq. 1) then
         nl(i) = n_str//l_str(ll+1)//' '
      else
         nl(i) = n_str//l_str(ll+1)//'-'
      end if

      ! check grid
      if (npts .gt. npts0) then
         write (*, *) 'error: npts .gt. npts0'
         stop
      end if

      ! read radial functions and grid
      read (3) a0, (pg(j,i), j=1, npts), (qg(j,i), j=1, npts)
      read (3) (rg(j,i), j=1, npts)

      ! find radial grid of longest extent and use that as a common grid
      ! in the output file
      if (i .gt. 1) then
         if ( maxval(rg(:,i)) .gt. maxval(rg(:,i-1)) ) then
            i_rg_max = i
         end if
      end if

      write(*, '(i5, a7, es14.6, i7, f10.2, es16.6)') i, nl(i), energy, npts, maxval(rg(:,i)), a0

      ! increase orbital counter
      i = i + 1

   end do

20 continue

   ! number of orbitals
   nwf = i - 1

   ! write header
   write (4, '(a14)', advance='no') 'r(a.u)'
   do i = 1, nwf
         !if ( any(selected_orbs == trim(adjustl(nl(i)))) ) then
         write (4, '(a14)', advance='no') 'P('//trim(adjustl(nl(i)))//')'
         write (4, '(a14)', advance='no') 'Q('//trim(adjustl(nl(i)))//')'
         !end if
   end do
   write (4, *)

   ! write r, P(r) and Q(r) for all orbitals
   do j = 1, npts
      write (4, '(es14.6)', advance='no') rg(j, nwf)
      do i = 1, nwf
         !if ( any(selected_orbs == trim(adjustl(nl(i)))) ) then
         write (4, '(2es14.6)', advance='no') pg(j,i), qg(j,i)
         !end if
      end do
      write (4, *)
   end do

   write(*,*)
   write(*,*) ' Done! Radial wavefunctions printed to rwfn.plot.'
   write(*,*)

end program
