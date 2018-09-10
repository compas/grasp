!     last edited August 2, 1996
      subroutine matbin(org, lock, closed, varmax, skal, second, anel0, par0, &
         low, nmax, lim, dubbel, minj, maxj) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax 
      integer , intent(out) :: skal 
      integer , intent(in) :: anel0 
      integer , intent(out) :: par0 
      integer  :: nmax 
      integer  :: minj 
      integer  :: maxj 
      logical , intent(inout) :: second 
      integer  :: org(15,0:10) 
      integer , intent(inout) :: low(15,0:10) 
      integer  :: lim(15) 
      logical , intent(inout) :: lock(15,0:10) 
      logical , intent(inout) :: closed(15,0:10) 
      logical , intent(out) :: dubbel(15,0:10) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: logfil = 31 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: anel, par, i, j, lmax, em, nenter, anela, anelb, block, mshell&
         , tmp 
      logical :: all, lima 
      character :: x 
      character , dimension(0:10) :: orb 
      character , dimension(0:20) :: l 
      character :: y*2 
!-----------------------------------------------
!
      data (l(i),i=0,20)/ 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N'&
         , 'O', 'Q', 'R', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'/  
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
 
   40 continue 
      if (.not.second) then 
         write (*, 200) 'Generate another list? (y/*)' 
         read (193, 1000) x 
         second = x=='y' .or. x=='Y' 
         write (logfil, *) second, ' Generate another list.' 
         if (.not.second) return  
      endif 
      anel = 0 
      anela = 0 
      anelb = 0 
      par = 0 
      skal = 20 
   60 continue 
      write (*, 200) 'Highest n-number? (1..15)' 
      read (193, *, err=60) nmax 
      nmax = max(nmax,1) 
      nmax = min(nmax,15) 
      write (logfil, *) nmax, ' Highest principal quantum number.' 
   70 continue 
      write (*, 400) 'Highest l-number? (s..', orb(min(10,nmax-1)), ')' 
      read (193, 1000) x 
      lmax = -1 
      do i = 0, min(10,nmax - 1) 
         if (x /= orb(i)) cycle  
         lmax = i 
      end do 
      if (lmax == (-1)) go to 70 
      write (logfil, *) lmax, ' Highest orbital angular momentum.' 
      write (*, 200) 'Are all these nl-subshells active? (n/*)' 
      read (193, 1000) x 
      all = .not.(x=='n' .or. x=='N') 
      write (logfil, *) all, ' all subshells active.' 
      lim = 0 
      if (nmax >= 2) then 
         write (*, 200) 'Limitations on population of n-subshells? (y/*)' 
         read (193, 1000) x 
         lima = x=='y' .or. x=='Y' 
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
               read (193, *, err=83) lim(i) 
               lim(i) = min0(mshell,lim(i)) 
               write (logfil, *) lim(i), &
                  ' is minimum number of electrons with n =', i 
            end do 
         endif 
      endif 
   95 continue 
      if (nmax < 10) then 
         write (*, 200) 'Highest n-number in reference configuration? (1..', &
            nmax, ')' 
      else 
         write (*, 202) 'Highest n-number in reference configuration? (1..', &
            nmax, ')' 
      endif 
      read (193, *, err=95) nenter 
      nenter = max(nenter,1) 
      nenter = min(nenter,nmax) 
      write (logfil, *) nenter, ' highest n-number.' 
      block = 0 
      do i = 1, 15 
         do j = 0, min(10,i - 1) 
            low(i,j) = 0 
            dubbel(i,j) = .FALSE. 
            if (nmax>=i .and. lmax>=j .and. .not.closed(i,j)) then 
               if (nenter >= i) then 
                  em = 2 + 4*j 
                  if (em < 10) then 
  100                continue 
                     if (i <= 9) then 
                        write (*, 200) 'Number of electrons in ', i, orb(j), &
                           '? (0..', em, ')' 
                     else 
                        write (*, 202) 'Number of electrons in ', i, orb(j), &
                           '? (0..', em, ')' 
                     endif 
                     read (193, *, err=100) org(i,j) 
                     if (org(i,j)<0 .or. org(i,j)>em) go to 100 
                  else 
  101                continue 
                     if (i < 10) then 
                        write (*, 201) 'Number of electrons in ', i, orb(j), &
                           '? (0..', em, ')' 
                     else 
                        write (*, 203) 'Number of electrons in ', i, orb(j), &
                           '? (0..', em, ')' 
                     endif 
                     read (193, *, err=101) org(i,j) 
                     if (org(i,j)<0 .or. org(i,j)>em) go to 101 
                  endif 
                  write (logfil, *) org(i,j), ' number of electrons in', i, orb&
                     (j) 
                  if (all) then 
                     lock(i,j) = .FALSE. 
                  else 
                     if (org(i,j) > 1) then 
                        if (org(i,j) <= 10) then 
                           write (*, 201) &
                              'Inactive, active or minimum? (i/*/0..', org(i,j)&
                               - 1, ')' 
                        else 
                           write (*, 202) &
                              'Inactive, active or minimum? (i/*/0..', org(i,j)&
                               - 1, ')' 
                        endif 
                        read (193, 1000) y 
                     else if (org(i,j) == 1) then 
                        write (*, 201) 'Inactive or active? (i/*)' 
                        read (193, 1000) y 
                     else 
                        write (*, 400) 'Inactive, active or doubled  ', &
                           'excited? (i/*/d)' 
                        read (193, 1000) y 
                        dubbel(i,j) = y(1:1)=='d' .or. y(1:1)=='D' 
                     endif 
                     lock(i,j) = y(1:1)=='i' .or. y(1:1)=='I' 
                     if (y(1:1)>='0' .and. y(1:1)<='9') then 
                        if (org(i,j) > 0) then 
                           tmp = ichar(y(1:1)) - ichar('0') 
                           if (y(2:2)>='1' .and. y(2:2)<='9') tmp = tmp*10 + &
                              ichar(y(2:2)) - ichar('0') 
                           low(i,j) = min(org(i,j),tmp) 
                        endif 
                     endif 
                     write (logfil, 1000) y, ' inactive, active, etc...' 
                  endif 
                  if (.not.lock(i,j)) anela = anela + org(i,j) 
                  anel = anel + org(i,j) 
                  par = mod(par + j*org(i,j),2) 
               else if (all) then 
                  org(i,j) = 0 
                  lock(i,j) = .FALSE. 
               else 
                  org(i,j) = 0 
                  closed(i,j) = .FALSE. 
                  if (i < 10) then 
                     write (*, 204) i, orb(j), ' inactive, active or ', &
                        'doubled excited? (i/*/d)' 
                  else 
                     write (*, 205) i, orb(j), ' inactive, active or ', &
                        'doubled excited? (i/*/d)' 
                  endif 
                  read (193, 1000) x 
                  dubbel(i,j) = x=='d' .or. x=='D' 
                  lock(i,j) = x=='i' .or. x=='I' 
                  write (logfil, *) x, i, orb(j), ' inactive, active, etc...' 
               endif 
            else 
               org(i,j) = 0 
               lock(i,j) = .TRUE. 
               if (closed(i,j)) then 
                  if (i < 10) then 
                     write (*, 204) i, orb(j), ' is a closed shell.' 
                  else 
                     write (*, 205) i, orb(j), ' is a closed shell.' 
                  endif 
                  em = 2 + 4*j 
                  anel = anel + em 
                  block = block + em 
               endif 
            endif 
            anelb = anelb + low(i,j) 
         end do 
         lim(i) = lim(i) - block 
         lim(i) = max0(0,lim(i)) 
      end do 
      if (anel /= anel0) then 
         if (anel0 < 10) then 
            write (*, 300) 'Wrong number of electrons. The first list had ', &
               anel0, ' electrons.' 
         else 
            write (*, 301) 'Wrong number of electrons. The first list had ', &
               anel0, ' electrons.' 
         endif 
         if (anel < 10) then 
            write (*, 300) 'This list has ', anel, ' electrons.' 
         else 
            write (*, 301) 'This list has ', anel, ' electrons.' 
         endif 
         second = .FALSE. 
         go to 40 
      endif 
 1100 continue 
      write (*, 400) 'Resulting 2*J-number? lower, higher ', &
         '(J=1 -> 2*J=2 etc.)' 
      read (193, *, err=1100) minj, maxj 
      if (anel == 2*(anel/2)) then 
         if (minj/=2*(minj/2) .or. maxj/=2*(maxj/2)) then 
            write (*, *) 'The resulting 2*J-numbers should be even' 
            go to 1100 
         endif 
      else 
         if (minj==2*(minj/2) .or. maxj==2*(maxj/2)) then 
            write (*, *) 'The resulting 2*J-numbers should be odd' 
            go to 1100 
         endif 
      endif 
      write (logfil, *) minj, ' to', maxj, ' is the resulting term.' 
!     if (par.NE.par0) then
!        write(*,200) 'Wrong parity.'
!        if (par0.EQ.0) write(*,*)
!    :           'The first list had even parity and this list has odd.'
!        if (par0.EQ.1) write(*,*)
!    :           'The first list had odd parity and this list has even.'
!        second = .FALSE.
!        goto 40
!     endif
      par0 = par 
      anelb = anela - anelb 
 1200 continue 
      if (anelb < 10) then 
         write (*, 200) 'Number of excitations = ? (0..', anelb, ')' 
         read (193, *, err=1200) varmax 
      else 
         write (*, 202) 'Number of excitations = ? (0..', anelb, ')' 
         read (193, *, err=1200) varmax 
      endif 
      write (logfil, *) varmax, ' number of excitations.' 
  200 format(' ',a,i1,a,a,i1,a) 
  201 format(' ',a,i1,a,a,i2,a) 
  202 format(' ',a,i2,a,a,i1,a) 
  203 format(' ',a,i2,a,a,i2,a) 
  204 format(' ',i1,3a) 
  205 format(' ',i2,3a) 
  208 format(' ',a,i1,a,i2,a) 
  300 format(' ',a,i1,a) 
  301 format(' ',a,i2,a) 
  400 format(' ',3a) 
  401 format(' ',2a,i1,a) 
  402 format(' ',2a,i2,a) 
 1000 format(a,a,a) 
 2000 format(i1,a) 
      return  
      end subroutine matbin 
