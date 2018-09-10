!***********************************************************************
!     ------------------------------------------------------------------
!     JJGEN  -- program to generate
!
!     Written by: Lennart Sturesson
!
!     2:nd version
!
!     Last edited Januar 2, 1997
!
!     ------------------------------------------------------------------
!
      program jjgen15 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reffa_I 
      use adder_I 
      use matcin_I 
      use blandc_I 
      use matain_I 
      use fivelines_I 
      use blanda_I 
      use matbin_I 
      use merge_I 
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: logfil = 31 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10) :: org 
      integer :: varmax, skal, anel, par 
      integer , dimension(15,0:10) :: low 
      integer , dimension(110) :: posn, posl 
      integer :: nmax 
      integer , dimension(15) :: lim 
      integer :: minj, maxj, cfmax 
      logical , dimension(15,0:10) :: lock, closed 
      logical :: slut 
      logical , dimension(15,0:10) :: med, dubbel 
      logical :: advexp, second 
      character :: x 
!-----------------------------------------------
      open(unit=logfil, file='clist.log', status='unknown', position='asis') 
      write (*, *) 'Version 2' 
      write (*, 200) '    *  : new list' 
      write (*, 200) '    a  : add to existing list' 
      write (*, 200) '    e  : expand existing list' 
      write (*, 200) '    q  : quit' 
      read (*, 100) x 
      write (logfil, 200) ' Option : ', x 
      call reffa (posn, posl) 
      advexp = .FALSE. 
      if (x=='a' .or. x=='A') then 
         call adder (closed, med, slut, anel, par, .FALSE.) 
         if (slut) then 
            write (*, 200) 
            write (*, 200) 'The clist.inp-file is not readable! ' 
            stop  
         endif 
         write (logfil, 200) ' New reference set.' 
         second = .TRUE. 
      else if (x=='e' .or. x=='E') then 
         advexp = .TRUE. 
         call adder (closed, med, slut, anel, par, .TRUE.) 
         if (slut) then 
            write (*, 200) 
            write (*, 200) 'The clist.inp-file is not readable! ' 
            stop  
         endif 
         write (logfil, 200) ' File as reference sets.' 
         call matcin (lock, closed, med, varmax, cfmax, nmax, minj, maxj, lim) 
         call blandc (varmax, cfmax, lock, med, minj, maxj, nmax, posn, posl, &
            lim) 
         second = .FALSE. 
      else 
         call matain (org, lock, closed, varmax, skal, nmax, anel, par, low, &
            minj, maxj, lim, dubbel) 
         call fivelines (org, lock, closed, .TRUE., posn, posl) 
         call blanda (org, varmax, lock, minj, maxj, skal, nmax, low, posn, &
            posl, lim, dubbel, .TRUE.) 
         second = .FALSE. 
      endif 
      if (.not.advexp) call matbin (org, lock, closed, varmax, skal, second, &
         anel, par, low, nmax, lim, dubbel, minj, maxj) 
      if (second) then 
         call fivelines (org, lock, closed, .FALSE., posn, posl) 
         call blanda (org, varmax, lock, minj, maxj, skal, nmax, low, posn, &
            posl, lim, dubbel, .FALSE.) 
         call merge (.FALSE., posn, posl) 
         write (*, 200) 'The merged file is called clist.out.' 
      else 
         call merge (.TRUE., posn, posl) 
         write (*, 200) 'The generated file is called clist.out.' 
      endif 
      stop  
  100 format(a) 
  200 format(' ',10a) 
  300 format(' ',a,i2,a) 
  400 format(' ',a,i3,a) 
      stop  
      end program jjgen15 
