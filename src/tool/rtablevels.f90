program renergytable

! This program makes ASCII and LaTeX tables over energies as a function of
! increasing active sets

! Give the energy lists from rlevels. Start with list corresponding 
! to the smallest active set and end with the one corresponding to the
! largest. Labels and order of energy levels are according to the last
! energy list. A maximum of 1000 levels in 15 files are allowed.

! Per Jonssson, Malmo University, August 2014

implicit none
integer :: h,i,j,k,l, nfile, nlevels, maxlengthascii, maxlengthlatex, nlength, ncount, ncase, nskip, nsame, lastpos
character(len=100) :: filename
character(len=200) :: fileline
character(len=1) :: extra(15,1000)
character(len=2) :: labelchar(1000)
character(len=8) :: jparity(15,1000)
character(len=12) :: energy(15,1000)
!character(len=12) :: energy(12,1000)
character(len=145) :: label(15,1000), labelstring, dummystring, latexstring(1000)
character(len=300) :: energystring, filenames
character(len=29) :: string1, string2, string3
character(len=1) :: char1, char2, char3
character(len=15) :: rstring
character(len=1) :: ans

open(unit=19,file='energytablelatex.tex',status='unknown')
open(unit=20,file='energytableascii.txt',status='unknown')

write(*,*)
write(*,*) ' RTABLEVELS'
write(*,*) ' Makes LaTeX and ASCII tables of energy files produced by' 
write(*,*) ' rlevels (in ljs format)                       '
write(*,*) ' Multiple energy files can be used as input'
write(*,*) ' Energies from file 1 fills column 1, energies from file 2'
write(*,*) ' fills column 2 etc.  Checks are done to see if the labels'
write(*,*) ' if the labels in the files are consistent'
write(*,*) ' Input file: name1, name2, ...'
write(*,*) ' Output files: energylabellatex.tex, energylabelascii.txt' 
write(*,*)
write(*,*) ' Inspect energy files and determine how many positions'
write(*,*) ' should be skipped in the string that determines the label '
write(*,*) ' e.g. if the string is 1s(2).2s_2S.2p(2)3P2_4P and 1s(2) is a core'
write(*,*) ' then you would like to skip 1s(2). i.e. 6 positions and determine'
write(*,*) ' the label from 2s_2S.2p(2)3P2_4P'
write(*,*)
write(*,*) ' How many positions should be skipped?'
read(*,*) nskip

write(*,*) ' Give the number of energy files from rlevels'
read(*,*) nfile

rstring = ''
do i = 1,nfile
   rstring = trim(rstring)//'r'
end do

write(19,'(a)') '\documentclass[10pt]{article}'
write(19,'(a)') '\usepackage{longtable}'
write(19,'(a)') '\begin{document}'
write(19,'(a)') '\begin{longtable}{l'//trim(rstring)//'} \hline'

jparity(:,:) = ''
labelchar(:) = ''
filenames = ''

nsame = 0
do i = 1,nfile
   write(*,'(a,i2)') '  Name of file',i
   read(*,'(a)') filename

   filenames = trim(filenames)//trim(filename)//', '

   open(unit=20+i,file=trim(filename),status='old')

! Start reading the file 

   k = 0
   do 
      read(20+i,'(a)') fileline
      if (fileline(1:3).eq.'---') then
         k = k + 1
         if (k.eq.2) exit
      end if
   end do

   nlevels = 0
   do
      read(20+i,'(a)') fileline
      if (fileline(1:3).eq.'---') exit
      nlevels = nlevels + 1
      jparity(i,nlevels) = fileline(7:14)
!      energy(i,nlevels) = fileline(30:41)
      energy(i,nlevels) = fileline(30:38)
!      label(i,nlevels) = fileline(56:200)
      label(i,nlevels) = fileline(56+nskip:200)
      lastpos = len_trim(label(i,nlevels))

      select case (label(i,nlevels)(lastpos:lastpos))
      case ('S','P','D','F','G','H','I','K','L','M','N')
         extra(i,nlevels) = ' ' 
      case default
         extra(i,nlevels) = label(i,nlevels)(lastpos:lastpos)
          label(i,nlevels)(lastpos:lastpos) = ' '
      end select
!      write(*,*) 'Label, jparity',trim(label(i,nlevels)),trim(jparity(i,nlevels)),extra(i,nlevels)
   end do 

! Check if any levels have the same J, parity and label

   do j = 1,nlevels
      do k = j+1,nlevels
         if (jparity(i,j).eq.jparity(i,k).and.label(i,j).eq.label(i,k)) then
            nsame = nsame + 1
            write(*,*) 'Level',j,'and level',k,'in file',i,'have the same label'
         end if
      end do
   end do

end do

if (nsame.gt.0) then
   write(*,*) 
   write(*,*)
   write(*,*) 'There are levels with the same labels extra character added'
   write(*,*) 'at the end to get unique labels'
   
   i = nfile
   do j = 1,nlevels
      nsame = 0
      do k = j+1,nlevels
         if (jparity(i,j).eq.jparity(i,k).and.label(i,j).eq.label(i,k).and.labelchar(k).eq.'  ') then
            nsame = nsame + 1
            if (nsame.eq.1) then
               labelchar(k) = '~b'
            elseif (nsame.eq.2) then
               labelchar(k) = '~c'
            elseif (nsame.eq.3) then
               labelchar(k) = '~d'
            else
               write(*,*) 'Too many labels are equal'
               stop
            end if
         end if
      end do
      if (nsame.gt.0) labelchar(j) = '~a'
   end do
end if

! Convert quantum labels for the last file (corresponding to the largest active set)
! to LaTeX

do h = 1,nlevels
   labelstring = label(nfile,h)

   do i = 1,142

!  Replace (n) with ^n

      if ((labelstring(i:i).eq.'(').and.(labelstring(i+2:i+2).eq.')')) then
         labelstring(i:i) = '^'
         labelstring(i+2:i+2) = ' '
      end if
   end do

   do i = 1,142

!  Replace . with \,

      if (labelstring(i:i).eq.'.') then
         dummystring = labelstring
         labelstring(1:i-1) = dummystring(1:i-1)
         labelstring(i:i) = '\'
         labelstring(i+1:i+1) = ','
         labelstring(i+2:145) = dummystring(i+1:143) 
      end if
   end do

   do i = 1,142

!  Replace _ with ~

      if (labelstring(i:i).eq.'_') labelstring(i:i) = '~'
   end do
!   write(*,'(a)') trim(labelstring)

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and integer2 replace with (^integer1_integer2S), (^integer1_integer2P), etc

   do l = 1,15
      ncase = 0
      do i = 1,142
         do j = 48,57
            do k = 48,57 
               char1 = labelstring(i:i)
               char2 = labelstring(i+1:i+1)
               char3 = labelstring(i+2:i+2)
               if ((ichar(char1).eq.j).and.(ichar(char3).eq.k).and.((char2.ne.'~').and.(char2.ne.' ').and.(char2.ne.'_'))) then
                  dummystring = labelstring
                  labelstring(1:i-1) = dummystring(1:i-1)
                  labelstring(i:i+6) = '(^'//char1//'_'//char3//char2//')'
                  labelstring(i+7:145) = dummystring(i+3:141)
                  ncase = ncase + 1
               end if
            end do
         end do
         if (ncase.eq.1) exit
      end do

!      write(*,'(a)') trim(labelstring)
   end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and not integer2 replace with ^integer1S, ^integer1P, etc 

   do i = 1,142
!  
      if (labelstring(i:i).eq.'~') then
         dummystring = labelstring
         labelstring(1:i) = dummystring(1:i)
         labelstring(i+1:i+1) = '^'
         labelstring(i+2:145) = dummystring(i+1:143)
      end if
   end do

   if  (jparity(nfile,h)(7:7).eq.'-') then 
      latexstring(h) = '$'//trim(labelstring)//'_{'//jparity(nfile,h)(1:5)//labelchar(h)//'}^o$'
   else
      latexstring(h) = '$'//trim(labelstring)//'_{'//jparity(nfile,h)(1:5)//labelchar(h)//'}$'
   end if
   
!  write(19,'(a)') trim(latexstring(h))

end do


! Start matching labels together for the different orbital sets

maxlengthascii = 0
maxlengthlatex = 0

do i = 1,nlevels
   nlength = len_trim(label(nfile,i))
   if (nlength.gt.maxlengthascii) maxlengthascii = nlength
   nlength = len_trim(latexstring(i))
   if (nlength.gt.maxlengthlatex) maxlengthlatex = nlength
end do

do i = 1,nlevels
   energystring = ''
   do j = 1,nfile-1
      ncount = 0
      do k = 1,nlevels
         if ((jparity(j,k).eq.jparity(nfile,i)).and.(label(j,k).eq.label(nfile,i))) then
            ncount = ncount + 1
            if (ncount.eq.1) then
              energystring = trim(energystring)//' & '//energy(j,k)
            end if
         end if
      end do
      if (ncount.eq.0) then
         write(*,*) 'Labels of levels are those of the last file'
         write(*,*) 'There are problems to merge data from all files to a table'
         write(*,*) 'Level',i,'in the last file is not found in file',j
         stop
      elseif (ncount.gt.1) then
         write(*,*) 'Labels of levels are those of the last file'
         write(*,*) 'Level',i,'in the last file occurs more than once in file',j
         write(*,*) 'Check that the merging is correct'
      end if
   end do
!   write(19,'(a)') latexstring(i)(1:maxlengthlatex)//trim(energystring)//' & '//energy(nfile,i)//' \\'
   write(20,'(a)') jparity(nfile,i)//' '//label(nfile,i)(1:maxlengthascii)//extra(nfile,i)//trim(energystring)//' '//energy(nfile,i) 
!   write(20,'(a)') label(nfile,i)(1:maxlengthascii)//extra(nfile,i)//' '//jparity(nfile,i)//trim(energystring)//' '//energy(nfile,i) 
   write(19,'(a)') latexstring(i)(1:maxlengthlatex)//'~'//extra(nfile,i)//trim(energystring)//' & '//energy(nfile,i)//' \\'
end do

write(19,'(a)') '\hline\\'
write(19,'(a)') '\caption{Energies from the files '//trim(filenames)//'}'
write(19,'(a)') '\end{longtable}'
write(19,'(a)') '\end{document}'

end program renergytable


