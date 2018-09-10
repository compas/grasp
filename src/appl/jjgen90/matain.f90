!     last edited November 1, 1996
      subroutine matain(org, lock, closed, varmax, skal, nmax, anel, par, low, &
         minj, maxj, lim, dubbel) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:24:37   1/ 2/07  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax 
      integer , intent(out) :: skal 
      integer  :: nmax 
      integer , intent(out) :: anel 
      integer , intent(out) :: par 
      integer  :: minj 
      integer  :: maxj 
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
      integer :: anela, anelb, resl, i, j, lmax, em, nenter, block, mshell, enn&
         , tmp, cormax 
      logical :: log_all, lima, open_c, clos_c 
      character :: x 
      character , dimension(0:10) :: orb 
      character , dimension(0:20) :: l 
      character :: y*3 
!-----------------------------------------------
      data (l(i),i=0,20)/ 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N'&
         , 'O', 'Q', 'R', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'/  
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
 
      do i = 1, 15 
         org(i,:min(i-1,10)) = 0 
      end do 
      skal = 20 
   60 continue 
      write (*, 200) 'Highest principal quantum number, n? (1..15)' 
      read (*, *, err=60) nmax 
      nmax = max(nmax,1) 
      nmax = min(nmax,15) 
      write (logfil, *) nmax, ' Highest principal quantum number.' 
   70 continue 
      write (*, 300) 'Highest orbital angular momentum, l? (s..', orb(min(10,&
         nmax-1)), ')' 
      read (*, 1000) x 
      lmax = -1 
      do i = 0, min(10,nmax - 1) 
         if (x /= orb(i)) cycle  
         lmax = i 
      end do 
      if (lmax == (-1)) go to 70 
      write (logfil, *) lmax, ' Highest orbital angular momentum.' 
      write (*, 200) 'Are all these nl-subshells active? (n/*)' 
      read (*, 1000) x 
      log_all = .not.(x=='n' .or. x=='N') 
      write (logfil, *) log_all, ' all subshells active.' 
      lim = 0 
      if (nmax >= 2) then 
         write (*, 200) 'Limitations on population of n-subshells? (y/*)' 
         read (*, 1000) x 
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
               read (*, *, err=83) lim(i) 
               lim(i) = min0(mshell,lim(i)) 
               write (logfil, *) lim(i), &
                  ' is minimum number of electrons with n =', i 
            end do 
         endif 
      endif 
   90 continue 
      if (nmax < 10) then 
         write (*, 200) 'Highest n-number in reference configuration? (1..', &
            nmax, ')' 
      else 
         write (*, 202) 'Highest n-number in reference configuration? (1..', &
            nmax, ')' 
      endif 
      read (*, *, err=90) nenter 
      nenter = max(nenter,1) 
      nenter = min(nenter,nmax) 
      write (logfil, *) nenter, ' highest n-number.' 
      if (nenter >= 1) then 
         write (*, 200) 'Predefine open, closed or no core? (o/c/*)' 
         read (*, 1000) x 
         open_c = x=='O' .or. x=='o' 
         clos_c = x=='C' .or. x=='c' 
         write (logfil, *) 'Predefined core:', x 
         if (open_c .or. clos_c) then 
   92       continue 
            write (*, 200) 'Select core,' 
            write (*, 200) &
               '     1: He (       1s(2)                  =  2 electrons)' 
            if (nenter >= 2) write (*, 200) &
               '     2: Ne ([He] + 2s(2)2p(6)             = 10 electrons)' 
            if (nenter >= 3) write (*, 200) &
               '     3: Ar ([Ne] + 3s(2)3p(6)             = 18 electrons)' 
            if (nenter >= 4) write (*, 200) &
               '     4: Kr ([Ar] + 3d(10)4s(2)4p(6)       = 36 electrons)' 
            if (nenter >= 5) write (*, 200) &
               '     5: Xe ([Kr] + 4d(10)5s(2)5p(6)       = 54 electrons)' 
            if (nenter >= 6) write (*, 200) &
               '     6: Rn ([Xe] + 4f(16)5d(10)6s(2)6p(6) = 86 electrons)' 
            read (*, *, err=92) cormax 
            if (cormax > nenter) go to 92 
            if (cormax >= 1) then 
               write (logfil, *) 'Core n=', cormax 
               if (clos_c) then 
                  do i = 1, cormax 
                     do j = 0, min(3,i - 1) 
                        closed(i,j) = .TRUE. 
                        org(i,j) = 2 + 4*j 
                     end do 
                  end do 
               else 
                  do i = 1, cormax 
                     do j = 0, min(3,i - 1) 
                        org(i,j) = 2 + 4*j 
                     end do 
                  end do 
               endif 
               select case (cormax)  
               case (3)  
                  org(3,2) = 0 
                  if (clos_c) closed(3,2) = .FALSE. 
               case (4)  
                  org(4,2) = 0 
                  org(4,3) = 0 
                  if (clos_c) then 
                     closed(4,2) = .FALSE. 
                     closed(4,3) = .FALSE. 
                  endif 
               case (5)  
                  org(4,3) = 0 
                  org(5,2) = 0 
                  org(5,3) = 0 
                  if (clos_c) then 
                     closed(4,3) = .FALSE. 
                     closed(5,2) = .FALSE. 
                     closed(5,3) = .FALSE. 
                  endif 
               case (6)  
                  org(5,3) = 0 
                  org(6,2) = 0 
                  org(6,3) = 0 
                  if (clos_c) then 
                     closed(5,3) = .FALSE. 
                     closed(6,2) = .FALSE. 
                     closed(6,3) = .FALSE. 
                  endif 
               end select 
            else 
               write (logfil, *) 'Core cancelled' 
            endif 
         endif 
      endif 
      anela = 0 
      anelb = 0 
      anel = 0 
      par = 0 
      block = 0 
      do i = 1, 15 
         do j = 0, min(10,i - 1) 
            low(i,j) = 0 
            dubbel(i,j) = .FALSE. 
            if (nmax>=i .and. lmax>=j .and. org(i,j)==0) then 
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
                     read (*, *, err=100) org(i,j) 
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
                     read (*, *, err=101) org(i,j) 
                     if (org(i,j)<0 .or. org(i,j)>em) go to 101 
                  endif 
                  write (logfil, *) org(i,j), ' number of electrons in', i, orb&
                     (j) 
                  anel = anel + org(i,j) 
                  par = mod(par + j*org(i,j),2) 
                  if (log_all) then 
                     lock(i,j) = .FALSE. 
                     closed(i,j) = .FALSE. 
                  else 
                     if (org(i,j) == em) then 
                        if (org(i,j) <= 10) then 
                           write (*, 201) &
                              'Closed, inactive, active or minimum? (c/i/*/0..'&
                              , org(i,j) - 1, ')' 
                        else 
                           write (*, 202) &
                              'Closed, inactive, active or minimum? (c/i/*/0..'&
                              , org(i,j) - 1, ')' 
                        endif 
                        read (*, 1000) y 
                        write (logfil, *) y, ' closed, inactive, etc...' 
                        closed(i,j) = y(1:1)=='c' .or. y(1:1)=='C' 
                        lock(i,j) = y(1:1)=='i' .or. y(1:1)=='I' .or. closed(i,&
                           j) 
                        if (closed(i,j)) block = block + em 
                     else 
                        if (org(i,j) > 1) then 
                           if (org(i,j) <= 10) then 
                              write (*, 201) &
                                 'Inactive, active or minimum? (i/*/0..', org(i&
                                 ,j) - 1, ')' 
                           else 
                              write (*, 202) &
      'Closed, inactive, active or minimum? (c/i/*/0..', org(i,j) - 1, ')' 
                           endif 
                        else if (org(i,j) == 1) then 
                           write (*, 201) 'Inactive or active? (i/*)' 
                        else 
                           write (*, 400) 'Inactive, active or double ', &
                              'excited? (i/*/d)' 
                        endif 
                        read (*, 1000) y 
                        write (logfil, *) y, ' inactive, active, etc...' 
                        if (org(i,j) == 0) dubbel(i,j) = y(1:1)=='d' .or. y(1:1&
                           )=='D' 
                        lock(i,j) = y(1:1)=='i' .or. y(1:1)=='I' 
                        closed(i,j) = .FALSE. 
                     endif 
                     if (y(1:1)>='0' .and. y(1:1)<='9') then 
                        if (org(i,j) > 0) then 
                           tmp = ichar(y(1:1)) - ichar('0') 
                           if (y(2:2)>='1' .and. y(2:2)<='9') tmp = tmp*10 + &
                              ichar(y(2:2)) - ichar('0') 
                           low(i,j) = min(org(i,j),tmp) 
                        endif 
                     endif 
                  endif 
                  if (.not.lock(i,j)) anela = anela + org(i,j) 
               else if (log_all) then 
                  org(i,j) = 0 
                  lock(i,j) = .FALSE. 
                  closed(i,j) = .FALSE. 
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
                  read (*, 1000) x 
                  write (logfil, *) x, i, orb(j), ' inactive, active, etc...' 
                  dubbel(i,j) = x=='d' .or. x=='D' 
                  lock(i,j) = x=='i' .or. x=='I' 
               endif 
            else if (org(i,j) /= 0) then 
               write (*, 204) i, orb(j), ' is part of the predefined core.' 
               if (open_c) then 
                  if (log_all) then 
                     closed(i,j) = .FALSE. 
                     lock(i,j) = .FALSE. 
                  else 
                     if (org(i,j) <= 10) then 
                        write (*, 201) &
                           'Closed, inactive, active or minimum? (c/i/*/0..', &
                           org(i,j) - 1, ')' 
                     else 
                        write (*, 202) &
                           'Closed, inactive, active or minimum? (c/i/*/0..', &
                           org(i,j) - 1, ')' 
                     endif 
                     read (*, 1000) y 
                     write (logfil, *) y, ' closed, inactive, etc...' 
                     closed(i,j) = y(1:1)=='c' .or. y(1:1)=='C' 
                     lock(i,j) = y(1:1)=='i' .or. y(1:1)=='I' .or. closed(i,j) 
                     if (y(1:1)>='0' .and. y(1:1)<='9') then 
                        if (org(i,j) > 0) then 
                           tmp = ichar(y(1:1)) - ichar('0') 
                           if (y(2:2)>='1' .and. y(2:2)<='9') tmp = tmp*10 + &
                              ichar(y(2:2)) - ichar('0') 
                           low(i,j) = min(org(i,j),tmp) 
                        endif 
                     endif 
                  endif 
                  if (.not.lock(i,j)) anela = anela + org(i,j) 
               else 
                  lock(i,j) = closed(i,j) 
               endif 
               if (closed(i,j)) block = block + org(i,j) 
               anel = anel + org(i,j) 
            else 
               org(i,j) = 0 
               lock(i,j) = .TRUE. 
               closed(i,j) = .FALSE. 
            endif 
            anelb = anelb + low(i,j) 
         end do 
         lim(i) = lim(i) - block 
         lim(i) = max0(0,lim(i)) 
      end do 
 1100 continue 
      write (*, 400) 'Resulting 2*J-number? lower, higher ', &
         '(J=1 -> 2*J=2 etc.)' 
      read (*, *, err=1100) minj, maxj 
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
      anelb = anela - anelb 
 1200 continue 
      if (anelb < 10) then 
         write (*, 200) 'Number of excitations = ? (0..', anelb, ')' 
         read (*, *, err=1200) varmax 
      else 
         write (*, 202) 'Number of excitations = ? (0..', anelb, ')' 
         read (*, *, err=1200) varmax 
      endif 
      write (logfil, *) varmax, ' number of excitations.' 
  200 format(' ',a,i1,a,a,i1,a) 
  201 format(' ',a,i1,a,a,i2,a) 
  202 format(' ',a,i2,a,a,i1,a) 
  203 format(' ',a,i2,a,a,i2,a) 
  204 format(' ',i1,3a) 
  205 format(' ',i2,3a) 
  206 format(' ',i1,a,a,i2,a) 
  207 format(' ',i2,a,a,i2,a) 
  208 format(' ',a,i1,a,i2,a) 
  300 format(' ',3a) 
  400 format(' ',2a,i1,a) 
  402 format(' ',2a,i2,a) 
 1000 format(3a) 
 2000 format(i1,2a) 
 3000 format(a,i2,2a) 
      return  
      end subroutine matain 
