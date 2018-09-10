!     last edited August 1, 1996
      subroutine blandc(varmax, cfmax, lock, med, minj, maxj, nmax, posn, posl&
         , lim) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use blandb_I 
      use mergeb_I 
      use gen_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax 
      integer , intent(in) :: cfmax 
      integer  :: minj 
      integer  :: maxj 
      integer  :: nmax 
      integer  :: posn(110) 
      integer  :: posl(110) 
      integer  :: lim(15) 
      logical  :: lock(15,0:10) 
      logical , intent(in) :: med(15,0:10) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7 
      integer, parameter :: fil_2 = 8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: cf 
      integer , dimension(15,0:10) :: org 
      integer :: i, j, n, k, l, l1, antal, tal, antalc 
      integer , dimension(15,0:10) :: low 
      integer :: tot 
      integer , dimension(15000,15,0:10) :: lista 
      integer , dimension(15,0:10,0:1) :: ansats 
      integer :: par, start, stopp, skal, duplet, kvar 
      logical :: finns 
      logical, dimension(15000) :: lik 
      character :: rad*500 
      character, dimension(0:10) :: orb 
!-----------------------------------------------
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
         'n'/  
      cf = 0 
      antalc = 0 
      tot = 0 
      skal = 20 
      finns = .FALSE. 
      do i = 1, 15 
         org(i,:min(10,i-1)) = 0 
         low(i,:min(10,i-1)) = 0 
      end do 
      open(unit=7, status='scratch', position='asis') 
      read (fil_2, 1000) rad 
      write (fil_1, 1000) rad 
      read (fil_2, 1000) rad 
      write (fil_1, 1000) rad 
      read (fil_2, 1000) rad 
      write (fil_1, 1000) rad 
      read (fil_2, 1000) 
      do i = 1, 500 
         rad(i:i) = ' ' 
      end do 
      start = -2 
      stopp = 0 
      do k = 1, 110 
         i = posn(k) 
         j = posl(k) 
         if (.not.(med(i,j) .or. .not.lock(i,j))) cycle  
         start = start + 5 
         stopp = stopp + 5 
         rad(start:start) = char(ichar('0') + i) 
         rad(start+1:start+1) = orb(j) 
         if (j < 1) cycle  
         rad(start+2:start+2) = '-' 
         start = start + 5 
         stopp = stopp + 5 
         rad(start:start) = char(ichar('0') + i) 
         rad(start+1:start+1) = orb(j) 
      end do 
      write (fil_1, 1000) rad 
      read (fil_2, 1000) rad 
      write (fil_1, 1000) rad 
    2 continue 
      read (fil_2, 1000, end=200) rad 
      read (fil_2, 1000, end=200) 
      read (fil_2, 1000, end=200) 
      tot = tot + 1 
      do i = 1, 15 
         lista(tot,i,:min(10,i-1)) = 0 
      end do 
      do i = 1, skal 
         n = 9*(i - 1) + 3 
         l = n + 1 
         tal = ichar(rad(n:n)) - ichar('0') 
         if (tal>=1 .and. tal<=15) then 
            l1 = -1 
            do j = 0, tal - 1 
               if (orb(j) /= rad(l:l)) cycle  
               l1 = j 
            end do 
            if (l1 == (-1)) exit  
         else 
            exit  
         endif 
         antal = ichar(rad(l+3:l+3)) - ichar('0') 
         if (antal>=0 .and. antal<=9) then 
            antal = antal*10 
         else 
            antal = 0 
         endif 
         antal = antal + ichar(rad(l+4:l+4)) - ichar('0') 
         lista(tot,tal,l1) = lista(tot,tal,l1) + antal 
      end do 
      go to 2 
  200 continue 
      if (tot == 0) then 
         write (*, 1005) 'Nothing in inputfile!' 
         stop  
      endif 
      if (tot < 10) then 
         write (*, 2001) 'Number of csf:s in inputfile = ', tot 
      else if (tot < 100) then 
         write (*, 2002) 'Number of csf:s in inputfile = ', tot 
      else if (tot < 1000) then 
         write (*, 2003) 'Number of csf:s in inputfile = ', tot 
      else if (tot < 10000) then 
         write (*, 2004) 'Number of csf:s in inputfile = ', tot 
      else 
         write (*, 2005) 'Number of csf:s in inputfile = ', tot 
      endif 
      duplet = 0 
      lik(:tot) = .FALSE. 
      if (tot >= 2) then 
         do i = 1, tot - 1 
            if (lik(i)) cycle  
            l302: do j = i + 1, tot 
               if (lik(j)) cycle  l302 
               do k = 1, nmax 
                  do l = 0, min(10,k - 1) 
                     if (lista(i,k,l) == lista(j,k,l)) cycle  
                     cycle  l302 
                  end do 
               end do 
               lik(j) = .TRUE. 
            end do l302 
         end do 
      endif 
      duplet = duplet + count(lik(:tot)) 
      if (duplet < 10) then 
         write (*, 2001) 'Number of duplicat csf:s in file = ', duplet 
      else if (duplet < 100) then 
         write (*, 2002) 'Number of duplicat csf:s in file = ', duplet 
      else if (duplet < 1000) then 
         write (*, 2003) 'Number of duplicat csf:s in file = ', duplet 
      else if (duplet < 10000) then 
         write (*, 2004) 'Number of duplicat csf:s in file = ', duplet 
      else 
         write (*, 2005) 'Number of duplicat csf:s in file = ', duplet 
      endif 
      kvar = tot - duplet 
      write (*, *) 
      do i = 1, tot 
         if (.not.lik(i)) then 
            if (kvar > 1) then 
               if (kvar < 10) then 
                  write (*, 1001) kvar, ' csf:s still to be expanded.' 
               else if (kvar < 100) then 
                  write (*, 1002) kvar, ' csf:s still to be expanded.' 
               else if (kvar < 1000) then 
                  write (*, 1003) kvar, ' csf:s still to be expanded.' 
               else if (kvar < 10000) then 
                  write (*, 1004) kvar, ' csf:s still to be expanded.' 
               else 
                  write (*, 1006) kvar, ' csf:s still to be expanded.' 
               endif 
            else 
               write (*, 1005) 'The last csf is still to be expanded.' 
            endif 
            kvar = kvar - 1 
            par = 0 
            do k = 1, nmax 
               org(k,:min(k-1,10)) = lista(i,k,:min(k-1,10)) 
            end do 
            if (finns) then 
               open(unit=21, status='scratch', position='asis') 
               call blandb (org, nmax, varmax, lock, 21, low, lim, posn, posl, &
                  minj, maxj) 
               rewind (21) 
               call mergeb (antalc) 
               if (antalc < 10) then 
                  write (*, 2001) 'Number of uncoupled csf:s = ', antalc 
               else if (antalc < 100) then 
                  write (*, 2002) 'Number of uncoupled csf:s = ', antalc 
               else if (antalc < 1000) then 
                  write (*, 2003) 'Number of uncoupled csf:s = ', antalc 
               else if (antalc < 10000) then 
                  write (*, 2004) 'Number of uncoupled csf:s = ', antalc 
               else 
                  write (*, 2005) 'Number of uncoupled csf:s = ', antalc 
               endif 
            else 
               open(unit=20, status='scratch', position='asis') 
               call blandb (org, nmax, varmax, lock, 20, low, lim, posn, posl, &
                  minj, maxj) 
               rewind (20) 
               finns = .TRUE. 
               antalc = 0 
               write (*, 1005) 'The first configuration has been expanded.' 
            endif 
            if (antalc >= cfmax) then 
               write (*, 1005) 'Maximum number of uncoupled csf:s exceeded' 
               exit  
            endif 
         endif 
      end do 
      write (*, *) 
      write (*, 1005) 'Preparing the couplings of the csf:s.' 
 
      if (nmax < 15) then 
         do i = nmax + 1, 15 
            ansats(i,:min(10,i-1),0) = 0 
            ansats(i,:min(10,i-1),1) = 0 
         end do 
      endif 
      cf = 0 
  490 continue 
      do i = 1, 15 
         read (20, 5000, end=492) (ansats(i,j,0),j=0,min(10,i - 1)) 
         read (20, 5000, end=492) (ansats(i,j,1),j=0,min(10,i - 1)) 
      end do 
      par = 0 
      do i = 1, 15 
         do j = 0, min(10,i - 1) 
            do k = 0, min(j,1) 
               par = mod(par + j*ansats(i,j,k),2) 
            end do 
         end do 
      end do 
      call gen (ansats, posn, posl, skal, cf, .TRUE., minj, maxj, par) 
      go to 490 
  492 continue 
      rewind (fil_1) 
      if (cf == 0) then 
         write (*, 1005) 'No configuration state has been generated.' 
      else if (cf == 1) then 
         write (*, 1005) 'One configuration state has been generated.' 
      else if (cf < 10) then 
         write (*, 1001) cf, ' configuration states have been generated.' 
      else if (cf < 100) then 
         write (*, 1002) cf, ' configuration states have been generated.' 
      else if (cf < 1000) then 
         write (*, 1003) cf, ' configuration states have been generated.' 
      else if (cf < 10000) then 
         write (*, 1004) cf, ' configuration states have been generated.' 
      else if (cf < 100000) then 
         write (*, 1006) cf, ' configuration states have been generated.' 
      else 
         write (*, *) cf, ' configuration states have been generated.' 
      endif 
 1000 format(a) 
 1001 format(' ',i1,a) 
 1002 format(' ',i2,a) 
 1003 format(' ',i3,a) 
 1004 format(' ',i4,a) 
 1005 format(' ',a) 
 1006 format(' ',i5,a) 
 2001 format(' ',a,i1,'.') 
 2002 format(' ',a,i2,'.') 
 2003 format(' ',a,i3,'.') 
 2004 format(' ',a,i4,'.') 
 2005 format(' ',a,i5,'.') 
 
 5000 format(11i2) 
      return  
      end subroutine blandc 
