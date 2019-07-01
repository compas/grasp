program rasfsplit

! Written by Per Jonsson, Malmo University, November 2016

implicit none

integer :: npos,i,j,k,l,ios,nblock,nblockodd,nblockeven
integer :: nelec,ncftot,nw,nvectot,nvecsize
integer :: nb,nevblk(100),iatjp,iaspa,ivec(100)
integer :: system
integer, allocatable :: ncfblk(:)

double precision :: eav,eval(100)
double precision, allocatable :: evec(:,:)

character(len=100) :: name
character(len=700) :: string1, string2
character(len=7) :: blockstring(30),oddstring(15),evenstring(15)
character(len=6) :: G92MIX
character(len=1) :: ans

write(*,*) 'RASFSPLIT'
write(*,*) 'Splits an Atomic State Function made up by the files name.c, '
write(*,*) 'name.(c)m, name.w into the corresponding files for each '
write(*,*) 'parity and J block. If only the name.c file is available this'
write(*,*) 'file will be split'
write(*,*) 'Input files: name.c,name.(c)m, name.w'
write(*,*) 'Output files: name_even1.c, name_even1.(c)m. name_even1.w '
write(*,*) '              name_odd1.c, name_odd1.(c)m. name_odd1.w ...'
write(*,*)

write(*,*) 'Name of the state'
read(*,'(a)') name
write(*,*)

write(*,*) 'Each of the blocks must be built from the same orbital set'
write(*,*) 'This may not be true for MR expansions, but is normally true'
write(*,*) 'for SD-MR expansions'
write(*,*) 'Is the above condition fullfilled? (y,n)'
read(*,*) ans
write(*,*)
if ((ans.eq.'n').or.(ans.eq.'N')) then
   write(*,*) 'stop'
   stop
end if

! Initialize strings

oddstring(1)   = '_odd1  '
oddstring(2)   = '_odd2  '
oddstring(3)   = '_odd3  '
oddstring(4)   = '_odd4  '
oddstring(5)   = '_odd5  '
oddstring(6)   = '_odd6  '
oddstring(7)   = '_odd7  '
oddstring(8)   = '_odd8  '
oddstring(9)   = '_odd9  '
oddstring(10)  = '_odd10 '
oddstring(11)  = '_odd11 '
oddstring(12)  = '_odd12 '
oddstring(13)  = '_odd13 '
oddstring(14)  = '_odd14 '
oddstring(15)  = '_odd15 '

evenstring(1)  = '_even1 '
evenstring(2)  = '_even2 '
evenstring(3)  = '_even3 '
evenstring(4)  = '_even4 '
evenstring(5)  = '_even5 '
evenstring(6)  = '_even6 '
evenstring(7)  = '_even7 '
evenstring(8)  = '_even8 '
evenstring(9)  = '_even9 '
evenstring(10) = '_even10'
evenstring(11) = '_even11'
evenstring(12) = '_even12'
evenstring(13) = '_even13'
evenstring(14) = '_even14'
evenstring(15) = '_even15'

! Open name.c file

open(unit=36,file=trim(name)//'.c',status='old',iostat=ios)
if (ios.ne.0) then
   write(*,*) 'Problems with opening the file ',trim(name)//'.c'
   write(*,*) 'stop'
   stop
end if

! Read the file and determine which are the blocks

nblock     = 0
nblockodd  = 0
nblockeven = 0
do
   string1 = string2
   read(36,'(a)',end=99) string2
   if (string2(2:2).eq.'*') then
      nblock = nblock + 1
      npos = len_trim(string1)
      if (string1(npos:npos).eq.'-') then
         nblockodd = nblockodd + 1
         if (nblockodd.gt.15) then
            write(*,*) 'Number of odd blocks limited to 15'
            write(*,*) 'Modify the code for more'
            stop
         end if
         blockstring(nblock) = oddstring(nblockodd)
      else
         nblockeven = nblockeven + 1
         if (nblockeven.gt.15) then
            write(*,*) 'Number of even blocks limited to 15'
            write(*,*) 'Modify the code for more'
            stop
         end if
         blockstring(nblock) = evenstring(nblockeven)
      end if
   end if
end do

! Deal with the last block
99 continue

nblock = nblock + 1
npos = len_trim(string1)
if (string1(npos:npos).eq.'-') then
   nblockodd = nblockodd + 1
   blockstring(nblock) = oddstring(nblockodd)
else
   nblockeven = nblockeven + 1
   blockstring(nblock) = evenstring(nblockeven)
end if

write(*,*) 'nblock    ',nblock
write(*,*) 'nblockodd ',nblockodd
write(*,*) 'nblockeven',nblockeven
write(*,*)

! Now open the files and append the block designation

do i = 1,nblock
   open(unit=36+i,file=trim(name)//trim(blockstring(i))//'.c',status='unknown')
end do

! Rewind name.c, loop through and write to the different blocks

rewind (36)

do i = 1,5
   read(36,'(a)') string1
   do j = 1,nblock
      write(36+j,'(a)') trim(string1)
   end do
end do

i = 1
do
   read(36,'(a)',end=999) string1
   if (string1(2:2).eq.'*') then
      i = i + 1
   else
      write(36+i,'(a)') trim(string1)
   end if
end do

999 continue

close(36)

! Now open the name.(c)m files

do l = 1,2
   if (l.eq.1) then
      open(unit=36,file=trim(name)//'.m',status='old',form='unformatted',iostat=ios)
   else
      open(unit=36,file=trim(name)//'.cm',status='old',form='unformatted',iostat=ios)
   end if
   if (ios.eq.0) then
      if (l.eq.1) then
         write(*,*) 'File ',trim(name)//'.m  ','available'
      else
         write(*,*) 'File ',trim(name)//'.cm ','available'
      end if
      write(*,*)
      read(36,iostat=ios) G92MIX
      read(36) nelec, ncftot, nw, nvectot, nvecsize, nblock
      write(*,*) '  nelec    = ', nelec
      write(*,*) '  ncftot   = ', ncftot
      write(*,*) '  nw       = ', nw
      write(*,*) '  nvectot  = ', nvectot
      write(*,*) '  nvecsize = ', nvecsize
      write(*,*) '  nblock   = ', nblock
      write(*,*)

! Allocate needed arrays

      allocate( ncfblk(ncftot) )
      allocate( evec(100,ncftot) )

      write(*,*)
      write(*,*) 'Block data read from mixing file'
      write(*,*) '        block        ncf         nev        2j+1          parity'
      do i=1, nblock
         read(36,end=98) nb, ncfblk(i), nevblk(i), iatjp, iaspa
         write(*,*) nb, ncfblk(i), nevblk(i), iatjp, iaspa
         if (l.eq.1) then
            open(unit=36+i,file=trim(name)//trim(blockstring(i))//'.m',status='unknown',form='unformatted')
         else
            open(unit=36+i,file=trim(name)//trim(blockstring(i))//'.cm',status='unknown',form='unformatted')
         end if
! Write the first two lines
         write(36+i) G92MIX
         write(36+i) nelec, ncfblk(i), nw, nevblk(i), ncfblk(i)*nevblk(i), 1
! Write the data for the block
         write(36+i) 1, ncfblk(i), nevblk(i), iatjp, iaspa
         if(nevblk(i).gt.0) then
            read(36) (ivec(j), j = 1,nevblk(i))
            read(36) eav, (eval(j), j = 1, nevblk(i))
            read(36) ((evec(k,j), j = 1, ncfblk(i)), k = 1, nevblk(i))
            write(36+i) (ivec(j), j = 1,nevblk(i))
            write(36+i) eav, (eval(j), j = 1, nevblk(i))
            write(36+i) ((evec(k,j), j = 1, ncfblk(i)), k = 1, nevblk(i))
         end if
      end do
98 continue

      do i = 1,nblock
         close(36+i)
      end do

      deallocate(ncfblk)
      deallocate(evec)

   end if
   close(36)
end do

! Now we just copy name.w to the different blocks using shell script copy command

write(*,*)
do i = 1,nblock
   string1 = 'cp '//trim(name)//'.w  '//trim(name)//trim(blockstring(i))//'.w'
         j = system (trim(string1))
!cjb  use EXECUTE_COMMAND_LINE if function 'system' is not supported
!cjb  call execute_command_line (trim(string1), exitstat=j)
   write(*,*) 'Exit status of the name.w file copying for block',i,'was',j
end do

 end program rasfsplit
