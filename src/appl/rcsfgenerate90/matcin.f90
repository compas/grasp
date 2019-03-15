!     last edited July 31, 1996
      subroutine matcin(lock, closed, med, varmax, cfmax, nmax, minj, maxj, lim&
         )
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax
      integer  :: cfmax
      integer  :: nmax
      integer  :: minj
      integer  :: maxj
      integer  :: lim(15)
      logical , intent(out) :: lock(15,0:10)
      logical , intent(in) :: closed(15,0:10)
      logical , intent(in) :: med(15,0:10)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: logfil = 31
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, lmax, nmaks, mshell, lmaks, j
      logical :: all, lima
      character :: x
      character, dimension(0:10) :: orb
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/
      nmaks = 1
      lmaks = 0
      do i = 1, 15
         do j = 0, min(10,i - 1)
            if (.not.med(i,j)) cycle
            nmaks = i
            lmaks = max(j,lmaks)
         end do
      end do
   60 continue
      if (nmaks <= 9) then
         write (*, 201) 'Highest n-number? (', nmaks, '..15)'
      else
         write (*, 202) 'Highest n-number? (', nmaks, '..15)'
      endif
      read (*, *, err=60) nmax
      nmax = max(nmax,nmaks)
      nmax = min(nmax,15)
      write (logfil, *) nmax, ' Highest principal quantum number.'
      write (*, 400) 'Highest l-number? (', orb(lmaks), '..', orb(min(10,nmax-1&
         )), ')'
      read (*, 1000) x
      lmax = -1
      do i = 0, min(10,nmax - 1)
         if (x /= orb(i)) cycle
         lmax = i
      end do
      lmax = max(lmaks,lmax)
!     if (lmax.EQ.-1) goto 70
      write (logfil, *) lmax, ' Highest orbital angular momentum.'
      write (*, 200) 'Are all these nl-subshells active? (n/*)'
      read (*, 1000) x
      all = .not.(x=='n' .or. x=='N')
      write (logfil, *) all, ' all subshells active.'
      lim = 0
      if (nmax >= 2) then
!******************* modified by yu zou, 3/6/00
! this option cannot run correctly. It is not provided at present.
         lima = .FALSE.
!******************* modified by yu zou, 3/6/00
         write (logfil, *) lima, ' limitations on population of n-subshells.'
         if (lima) then
            mshell = 0
            do i = 1, nmax - 1
               mshell = mshell + 2*i*i
   83          continue
               if (i == 1) then
                  write (*, 200) 'Minimum number of electrons with n=1? (0..2)'
               else if (i < 10) then
                  if (mshell < 100) then
                     write (*, 208) 'Minimum number of electrons with n<=', i, &
                        '? (0..', mshell, ')'
                  else
                     write (*, 208) 'Minimum number of electrons with n<=', i, &
                        '? (0..)'
                  endif
               else
                  write (*, 202) 'Minimum number of electrons with n<=', i, &
                     '? (0..)'
               endif
               read (*, *, err=83) lim(i)
               lim(i) = min0(mshell,lim(i))
               write (logfil, *) lim(i), &
                  ' is minimum number of electrons with n =', i
            end do
         endif
      endif
      if (all) then
         do i = 1, 15
            if (i < 10) then
               do j = 0, min(10,i - 1)
                  if (nmax>=i .and. lmax>=j .and. .not.closed(i,j)) then
                     lock(i,j) = .FALSE.
                  else
                     lock(i,j) = .TRUE.
                     if (closed(i,j)) write (*, 204) i, orb(j), &
                        ' is a closed shell.'
                  endif
               end do
            else
               do j = 0, min(10,i - 1)
                  if (nmax>=i .and. lmax>=j .and. .not.closed(i,j)) then
                     lock(i,j) = .FALSE.
                  else
                     lock(i,j) = .TRUE.
                     if (closed(i,j)) write (*, 205) i, orb(j), &
                        ' is a closed shell.'
                  endif
               end do
            endif
         end do
      else
         do i = 1, 15
            if (i < 10) then
               do j = 0, min(10,i - 1)
                  if (nmax>=i .and. lmax>=j .and. .not.closed(i,j)) then
                     write (*, 204) i, orb(j), ' inactive or active? ', '(i/*)'
                     read (*, 1000) x
                     write (logfil, *) x, i, orb(j), &
                        ' inactive, active, etc...'
                     lock(i,j) = x=='i' .or. x=='I'
                  else
                     lock(i,j) = .TRUE.
                     if (closed(i,j)) write (*, 204) i, orb(j), &
                        ' is a closed shell.'
                  endif
               end do
            else
               do j = 0, min(10,i - 1)
                  if (nmax>=i .and. lmax>=j .and. .not.closed(i,j)) then
                     write (*, 205) i, orb(j), ' inactive or active? ', '(i/*)'
                     read (*, 1000) x
                     write (logfil, *) x, i, orb(j), &
                        ' inactive, active, etc...'
                     lock(i,j) = x=='i' .or. x=='I'
                  else
                     lock(i,j) = .TRUE.
                     if (closed(i,j)) write (*, 205) i, orb(j), &
                        ' is a closed shell.'
                  endif
               end do
            endif
         end do
      endif
 1100 continue
      write (*, 400) 'Resulting 2*J-number? lower, higher ', &
         '(J=1 -> 2*J=2 etc.)'
      read (*, *, err=1100) minj, maxj
      write (logfil, *) minj, ' to', maxj, ' is the resulting term.'
  160 continue
      write (*, 200) 'Number of excitations = ? (0..)'
      read (*, *, err=160) varmax
      write (logfil, *) varmax, ' number of excitations.'
  170 continue
      write (*, 400) 'Maximum number of uncoupled configuration', &
         ' states? (0..)'
      read (*, *, err=170) cfmax
      write (logfil, *) cfmax, ' maximum number '
      write (*, *)

  200 format(' ',a,i1,a,a,i1,a)
  201 format(' ',a,i1,a,a,i2,a)
  202 format(' ',a,i2,a,a,i1,a)
  203 format(' ',a,i2,a,a,i2,a)
  204 format(' ',i1,3a)
  205 format(' ',i2,3a)
  208 format(' ',a,i1,a,i2,a)
  300 format(' ',a,i1,a)
  301 format(' ',a,i2,a)
  400 format(' ',7a)
  401 format(' ',2a,i1,a)
  402 format(' ',2a,i2,a)
 1000 format(a,a,a)
 2000 format(i1,a)
      return
      end subroutine matcin
