
!==============================================================================
!     last edited November 1, 1996
      subroutine matain(org, lock, closed, varmax, skal, nmax, anel, par, low, &
         minj, maxj, lim, dubbel)
!==============================================================================
!     Editted by Charlotte Froese Fischer,    October, 2017
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
      logical :: all, log_all, lima, open_c, clos_c
      character :: x
      character , dimension(0:10) :: orb
      character , dimension(0:20) :: l
      character :: y*3

      data (l(i),i=0,20)/ 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N'&
         , 'O', 'Q', 'R', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'/
      data (orb(i),i=0,10)/ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', &
                            'n'/
!-----------------------------------------------

      do 10 i=1,15
         do 10 j=0,min(i-1,10)
   10       org(i,j) = 0
      skal = 20
   60 write(*,200) 'Highest principal quantum number, n? (1..15)'
      read(193,*,err=60) nmax
      nmax = max(nmax,1)
      nmax = min(nmax,15)
      write(logfil,*) nmax,' Highest principal quantum number.'
   70 write(*,300) 'Highest orbital angular momentum, l? (s..',       &
                                              orb(min(10,nmax-1)),')'
      read(193,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      if (lmax.EQ.-1) goto 70
      write(logfil,*) lmax,' Highest orbital angular momentum.'
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(193,1000) X
      all   = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(logfil,*) all,' all subshells active.'
      do 72 i=1,15
   72    lim(i) = 0
      if (nmax.GE.2) then
        write(*,200) 'Limitations on population of n-subshells? (y/*)'
         read(193,1000) X
         lima = X.EQ.'y' .OR. X.EQ.'Y'
         write(logfil,*) lima,         &
                          ' limitations on population of n-subshells.'
         if (lima) then
            mshell = 0
            do 85 i=1,nmax-1
               mshell = mshell + 2*i*i
   83          continue
               if (i.EQ.1) then
                  write(*,200)   &
                          'Minimum number of electrons with n=1? (0..2)'
               elseif (i.LT.10) then
                  if (mshell.LT.100) then
                     write(*,208) 'Minimum number of electrons with n<='&
                             ,i,'? (0..',mshell,')'
                  else
                     write(*,208) 'Minimum number of electrons with n<='&
                             ,i,'? (0..)'
                  endif
               else
                  write(*,202) 'Minimum number of electrons with n<=',&
                                i, '? (0..)'
               endif
               read(193,*,err=83) lim(i)
               if (lim(i).GT.mshell) lim(i) = mshell
               write(logfil,*) lim(i),                                   &
                           ' is minimum number of electrons with n =',i
   85       continue
         endif
      endif
   90 continue
      if (nmax.LT.10) then
         write(*,200)                                                  &
           'Highest n-number in reference configuration? (1..', nmax,')'
      else
         write(*,202)                                                  &
           'Highest n-number in reference configuration? (1..', nmax,')'
      endif
      read(193,*,err=90) nenter
      nenter = max(nenter,1)
      nenter = min(nenter,nmax)
      write(logfil,*) nenter,' highest n-number.'
      if (nenter.GE.1) then
         write(*,200) 'Predefine open, closed or no core? (o/c/*)'
         read(193,1000) X
         open_c = X.EQ.'O' .OR. X.EQ.'o'
         clos_c = X.EQ.'C' .OR. X.EQ.'c'
         write(logfil,*) 'Predefined core:',X
         if (open_c .OR. clos_c) then
   92       write(*,200) 'Select core,'
            write(*,200)                  &
            '     1: He (       1s(2)                  =  2 electrons)'
            if (nenter.GE.2) write(*,200) &
            '     2: Ne ([He] + 2s(2)2p(6)             = 10 electrons)'
            if (nenter.GE.3) write(*,200) &
            '     3: Ar ([Ne] + 3s(2)3p(6)             = 18 electrons)'
            if (nenter.GE.4) write(*,200) &
            '     4: Kr ([Ar] + 3d(10)4s(2)4p(6)       = 36 electrons)'
            if (nenter.GE.5) write(*,200) &
            '     5: Xe ([Kr] + 4d(10)5s(2)5p(6)       = 54 electrons)'
            if (nenter.GE.6) write(*,200) &
            '     6: Rn ([Xe] + 4f(16)5d(10)6s(2)6p(6) = 86 electrons)'
            read(193,*,err=92) cormax
	    if (cormax.GT.nenter) goto 92
            if (cormax.GE.1) then
               write(logfil,*) 'Core n=',cormax
	       do 93 i=1,cormax
                  do 93 j= 0,min(3,i-1)
		     if (clos_c) closed(i,j) = .TRUE.
   93                org(i,j) = 2+4*j
	       if (cormax.EQ.3) then
		  org(3,2) = 0
                  if (clos_c) closed(3,2) = .FALSE.
               elseif (cormax.EQ.4) then
                  org(4,2) = 0
                  org(4,3) = 0
                  if (clos_c) then
		     closed(4,2) = .FALSE.
                     closed(4,3) = .FALSE.
		  endif
               elseif (cormax.EQ.5) then
                  org(4,3) = 0
                  org(5,2) = 0
                  org(5,3) = 0
                  if (clos_c) then
                     closed(4,3) = .FALSE.
		     closed(5,2) = .FALSE.
                     closed(5,3) = .FALSE.
                  endif
               elseif (cormax.EQ.6) then
                  org(5,3) = 0
                  org(6,2) = 0
                  org(6,3) = 0
                  if (clos_c) then
                     closed(5,3) = .FALSE.
                     closed(6,2) = .FALSE.
                     closed(6,3) = .FALSE.
                  endif
               endif
            else
               write(logfil,*) 'Core cancelled'
            endif
         endif
      endif
      anela  = 0
      anelb  = 0
      anel   = 0
      par    = 0
      block  = 0
      do 160 i=1,15
         do 150 j=0,min(10,i-1)
            low(i,j)    = 0
            dubbel(i,j) = .FALSE.
            if (nmax.GE.i .AND. lmax.GE.j .AND. org(i,j).EQ.0) then
               if (nenter.GE.i) then
                  em = 2 + 4*j
                  if (em.LT.10) then
  100                continue
                     if (i.LE.9) then
                        write(*,200) 'Number of electrons in ',i,      &
                                                orb(j),'? (0..',em,')'
                     else
                        write(*,202) 'Number of electrons in ',i,      &
                                                orb(j),'? (0..',em,')'
                     endif
                     read(193,*,err=100) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em)            &
                              goto 100
                  else
  101                continue
                     if (i.LT.10) then
                        write(*,201) 'Number of electrons in ',i,      &
                                                 orb(j),'? (0..',em,')'
                     else
                        write(*,203) 'Number of electrons in ',i,      &
                                                 orb(j),'? (0..',em,')'
                     endif
                     read(193,*,err=101) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em)            &
                              goto 101
                  endif
                  write(logfil,*) org(i,j),' number of electrons in',i,&
                              orb(j)
                  anel = anel + org(i,j)
                  par  = mod(par+j*org(i,j),2)
                  if (all) then
                     lock(i,j)   = .FALSE.
                     closed(i,j) = .FALSE.
                  else
                     if (org(i,j).EQ.em) then
                        if (org(i,j).LE.10) then
                           write(*,201)                                &
                    'Closed, inactive, active or minimum? (c/i/*/0..', &
                                     org(i,j)-1,')'
                     else
                           write(*,202)                                &
                    'Closed, inactive, active or minimum? (c/i/*/0..', &
                                     org(i,j)-1,')'
                     endif
                        read(193,1000) Y
                        write(logfil,*) Y,' closed, inactive, etc...'
                        closed(i,j) = Y(1:1).EQ.'c' .OR. Y(1:1).EQ.'C'
                        lock(i,j)   = Y(1:1).EQ.'i' .OR. Y(1:1).EQ.'I' &
                                               .OR. closed(i,j)
                        if (closed(i,j)) block = block + em
                     else
                        if (org(i,j).GT.1) then
                           if (org(i,j).LE.10) then
                              write(*,201)                             &
                'Inactive, active or minimum? (i/*/0..',org(i,j)-1,')'
                           else
                            write(*,202)                              &
                    'Closed, inactive, active or minimum? (c/i/*/0..',&
                                     org(i,j)-1,')'
                           endif
                        elseif (org(i,j).EQ.1) then
                           write(*,201) 'Inactive or active? (i/*)'
                        else
                           write(*,400) 'Inactive, active or double ',&
                                  'excited? (i/*/d)'
                        endif
                        read(193,1000) Y
                        write(logfil,*) Y,' inactive, active, etc...'
                        if (org(i,j).EQ.0) then
                           dubbel(i,j) = Y(1:1).EQ.'d' .OR.            &
                                         Y(1:1).EQ.'D'
                        endif
                        lock(i,j)   = Y(1:1).EQ.'i' .OR. Y(1:1).EQ.'I'
                        closed(i,j) = .FALSE.
                     endif
                     if (Y(1:1).GE.'0' .AND. Y(1:1).LE.'9') then
                        if (org(i,j).GT.0) then
                           tmp = ICHAR(Y(1:1))-ICHAR('0')
                           if (Y(2:2).GE.'0' .AND. Y(2:2).LE.'9')      &
                              tmp = tmp*10 + ICHAR(Y(2:2))-ICHAR('0')
                           low(i,j) = min(org(i,j),tmp)
                        endif
                     endif
                  endif
                  if (.NOT. lock(i,j)) anela = anela + org(i,j)
               elseif (all) then
                  org(i,j) = 0
                  lock(i,j) = .FALSE.
                  closed(i,j) = .FALSE.
               else
                  org(i,j) = 0
                  closed(i,j) = .FALSE.
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive, active or ',    &
                              'doubled excited? (i/*/d)'
                  else
                     write(*,205) i,orb(j),' inactive, active or ',    &
                              'doubled excited? (i/*/d)'
                  endif
                  read(193,1000) X
                  write(logfil,*) X,i,orb(j),' inactive, active, etc...'
                  dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I'
               endif
            elseif (org(i,j).NE.0) then
               write(*,204) i,orb(j),' is part of the predefined core.'
               if (open_c) then
                  if (all) then
                     closed(i,j) = .FALSE.
                     lock(i,j)   = .FALSE.
                  else
                     if (org(i,j).LE.10) then
                        write(*,201)                                   &
                    'Closed, inactive, active or minimum? (c/i/*/0..', &
                                     org(i,j)-1,')'
                     else
                        write(*,202)                                   &
                    'Closed, inactive, active or minimum? (c/i/*/0..', &
                                     org(i,j)-1,')'
                     endif
                     read(193,1000) Y
                     write(logfil,*) Y,' closed, inactive, etc...'
                     closed(i,j) = Y(1:1).EQ.'c' .OR. Y(1:1).EQ.'C'
                     lock(i,j)   = Y(1:1).EQ.'i' .OR. Y(1:1).EQ.'I'    &
                                               .OR. closed(i,j)
                     if (Y(1:1).GE.'0' .AND. Y(1:1).LE.'9') then
                        if (org(i,j).GT.0) then
                           tmp = ICHAR(Y(1:1))-ICHAR('0')
                           if (Y(2:2).GE.'0' .AND. Y(2:2).LE.'9')      &
                              tmp = tmp*10 + ICHAR(Y(2:2))-ICHAR('0')
                           low(i,j) = min(org(i,j),tmp)
                        endif
                     endif
                  endif
                  if (.NOT. lock(i,j)) anela = anela + org(i,j)
               else
                  lock(i,j) = closed(i,j)
               endif
               if (closed(i,j)) block = block + org(i,j)
               anel = anel + org(i,j)
            else
               org(i,j)  = 0
               lock(i,j) = .TRUE.
               closed(i,j) = .FALSE.
            endif
            anelb = anelb + low(i,j)
  150       continue
            lim(i) = lim(i) - block
            if (lim(i).LT.0) lim(i) = 0
  160 continue
 1100 write(*,400) 'Resulting 2*J-number? lower, higher ',             &
                    '(J=1 -> 2*J=2 etc.)'
      read(193,*,ERR=1100) minJ,maxJ
      if (anel .EQ. 2*(anel/2)) then
         if (minJ .NE. 2*(minJ/2) .OR. maxJ .NE. 2*(maxJ/2)) then
            write(*,*) 'The resulting 2*J-numbers should be even'
            goto 1100
         endif
      else
         if (minJ .EQ. 2*(minJ/2) .OR. maxJ .EQ. 2*(maxJ/2)) then
            write(*,*) 'The resulting 2*J-numbers should be odd'
            goto 1100
         endif
      endif
      write(logfil,*) minJ,' to',maxJ,' is the resulting term.'
      anelb = anela - anelb
 1200 continue
      if (anelb.LT.10) then
         write(*,200) 'Number of excitations = ? (0..',anelb,')'
         read(193,*,err=1200) varmax
      else
         write(*,202) 'Number of excitations = ? (0..',anelb,')'
         read(193,*,err=1200) varmax
      endif
      write(logfil,*) varmax,' number of excitations.'
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  206 format(' ',I1,A,A,I2,A)
  207 format(' ',I2,A,A,I2,A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',3A)
  400 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(3A)
 2000 format(I1,2A)
 3000 format(A,I2,2A)
      return
      end
