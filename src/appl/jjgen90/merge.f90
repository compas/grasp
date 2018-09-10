!     last edited July 31, 1996
      subroutine merge(single, posn, posl, ii) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lika_I 
      use fivefirst_I 
      use lasa1_I 
      use test_I 
      use lasa2_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ii 
      logical , intent(in) :: single 
      integer  :: posn(110) 
      integer  :: posl(110) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7 
      integer, parameter :: fil_2 = 8 
      integer, parameter :: utfil = 9 
      integer, parameter :: nyfil = 13 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10,0:1) :: pop1, pop2 
      integer :: skal1, skal2, i, j, k, cf, stopp1, stopp2 
      integer , dimension(15,0:10,0:1) :: popo 
      integer :: ii1 
      logical :: p1, p2, slut1, slut2 
      character :: rad11*200, rad12*200, rad21*200, rad22*200, rad31*200, rad32&
         *200 
!-----------------------------------------------
 
      ii1 = mod(ii,2) 
      if (ii1 == 0) then 
         open(unit=utfil, file='clist.out', status='unknown', position='asis') 
      else 
         open(unit=utfil, file='fil1.dat', status='unknown', position='asis') 
      endif 
      open(unit=nyfil, file='clist.new', status='unknown', position='asis') 
      slut1 = .FALSE. 
      slut2 = single 
      cf = 0 
      call fivefirst (slut1, slut2, posn, posl) 
      skal1 = 20 
      skal2 = 20 
      call lasa1 (fil_1, rad11, pop1, skal1, slut1) 
      call lasa1 (fil_2, rad12, pop2, skal2, slut2) 
   10 continue 
      if (.not.slut1 .and. .not.slut2) then 
         call test (p1, p2, pop1, pop2, 15) 
         if (p1) then 
            do i = 1, 15 
               popo(i,:min(10,i-1),:1) = pop1(i,:min(10,i-1),:1) 
            end do 
            stopp1 = max(1,9*skal1) 
            stopp2 = 9*skal1 + 2 
   30       continue 
            call lasa2 (fil_1, rad21, rad31, stopp1, slut1) 
            if (.not.slut1) then 
               write (utfil, 999) rad11(1:stopp1) 
               write (utfil, 999) rad21(1:stopp1) 
               write (utfil, 999) rad31(1:stopp2) 
               cf = cf + 1 
            endif 
            skal1 = 20 
            call lasa1 (fil_1, rad11, pop1, skal1, slut1) 
            if (.not.slut1) then 
               if (lika(popo,pop1)) go to 30 
            endif 
            if (p2) then 
   40          continue 
               call lasa2 (fil_2, rad22, rad32, stopp1, slut2) 
               skal2 = 20 
               call lasa1 (fil_2, rad12, pop2, skal2, slut2) 
               if (.not.slut2) then 
                  if (lika(popo,pop2)) go to 40 
               endif 
            endif 
            go to 10 
         else if (p2) then 
            do i = 1, 15 
               popo(i,:min(10,i-1),:1) = pop2(i,:min(10,i-1),:1) 
            end do 
            stopp1 = max(1,9*skal2) 
            stopp2 = 9*skal2 + 2 
   60       continue 
            call lasa2 (fil_2, rad22, rad32, stopp1, slut2) 
            if (.not.slut2) then 
               write (utfil, 999) rad12(1:stopp1) 
               write (utfil, 999) rad22(1:stopp1) 
               write (utfil, 999) rad32(1:stopp2) 
               write (nyfil, 999) rad12(1:stopp1) 
               write (nyfil, 999) rad22(1:stopp1) 
               write (nyfil, 999) rad32(1:stopp2) 
               cf = cf + 1 
            endif 
            skal2 = 20 
            call lasa1 (fil_2, rad12, pop2, skal2, slut2) 
            if (.not.slut2) then 
               if (lika(popo,pop2)) go to 60 
            endif 
            go to 10 
         else 
            write (*, *) 'fatal error' 
            stop  
         endif 
      else if (.not.slut1 .and. slut2) then 
   70    continue 
         stopp1 = max(1,9*skal1) 
         stopp2 = 9*skal1 + 2 
         call lasa2 (fil_1, rad21, rad31, stopp1, slut1) 
         if (.not.slut1) then 
            write (utfil, 999) rad11(1:stopp1) 
            write (utfil, 999) rad21(1:stopp1) 
            write (utfil, 999) rad31(1:stopp2) 
            cf = cf + 1 
         endif 
         skal1 = 20 
         call lasa1 (fil_1, rad11, pop1, skal1, slut1) 
         if (.not.slut1) go to 70 
      else if (slut1 .and. .not.slut2) then 
   80    continue 
         stopp1 = max(1,9*skal2) 
         stopp2 = 9*skal2 + 2 
         call lasa2 (fil_2, rad22, rad32, stopp1, slut2) 
         if (.not.slut2) then 
            write (utfil, 999) rad12(1:stopp1) 
            write (utfil, 999) rad22(1:stopp1) 
            write (utfil, 999) rad32(1:stopp2) 
            write (nyfil, 999) rad12(1:stopp1) 
            write (nyfil, 999) rad22(1:stopp1) 
            write (nyfil, 999) rad32(1:stopp2) 
            cf = cf + 1 
         endif 
         skal2 = 20 
         call lasa1 (fil_2, rad12, pop2, skal2, slut2) 
         if (.not.slut2) go to 80 
      endif 
      close(fil_1) 
      close(fil_2) 
      close(utfil) 
      close(nyfil) 
      if (cf == 0) then 
         write (*, 105) 'No configuration state in the final list.' 
      else if (cf == 1) then 
         write (*, 105) 'One configuration state in the final list.' 
      else if (cf < 10) then 
         write (*, 101) cf, ' configuration states in the final list.' 
      else if (cf < 100) then 
         write (*, 102) cf, ' configuration states in the final list.' 
      else if (cf < 1000) then 
         write (*, 103) cf, ' configuration states in the final list.' 
      else if (cf < 10000) then 
         write (*, 104) cf, ' configuration states in the final list.' 
      else if (cf < 100000) then 
         write (*, 106) cf, ' configuration states in the final list.' 
      else 
         write (*, *) cf, ' configuration states in the final list.' 
      endif 
      return  
  101 format(' ',i1,a) 
  102 format(' ',i2,a) 
  103 format(' ',i3,a) 
  104 format(' ',i4,a) 
  105 format(' ',a) 
  106 format(' ',i5,a) 
  999 format(a) 
      return  
      end subroutine merge 
