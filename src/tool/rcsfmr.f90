program rcsfmr

! Written by er Jönsson
! Malmö University, Sweden, May 2018

implicit none
integer :: i,j,n
double precision :: cutoff,weight1,weight2
character(len=24) :: filename
character(len=500) :: longstring, string1,string2,string_vec(10000)

write(*,*) 'RCSFMR'
write(*,*) 'This program reads the name.lsj.lbl file and extracts a set of MR'
write(*,*) 'configuartions that give rise to LSJ coupled CSFs with absolute weights '
write(*,*) 'larger than a specified cut-off'
write(*,*) 'Input file: namel.lsj.lbl'
write(*,*) 'Ouput is written to screen'
write(*,*)
write(*,*) 'Name of state'
read(*,'(a)') filename
write(*,*) 'Give cut-off for weight'
read(*,*) cutoff

open(12,file=trim(filename)//'.lsj.lbl',form='formatted',status='old')
read(12,'(a)') longstring

i = 0
do
   read(12,'(a)',end=99) longstring
   n = len(trim(longstring))
   if ((n.ne.0).and.(longstring(n:n).ne.'%')) then
!      write(*,'(a)') trim(longstring)
      backspace 12
      read(12,*) weight1,weight2,string1
!      write(*,*) weight1,weight2
      if (dabs(weight1).gt.cutoff) then
         call stringprocess(string1,string2)
         i = i + 1
         string_vec(i) = string2
!         write(196,'(i5,a)') i,string_vec(i)
      end if
   end if
end do

99 continue

n = i ! Number of vectors

do i = 1,n
   do j = i+1,n
      if (trim(string_vec(i)).eq.trim(string_vec(j))) then
         string_vec(j) = ' '
      end if
   end do
end do

write(*,*)
write(*,*) 'Configurations in the MR'
write(*,*)

do i = 1,n
   if (len(trim(string_vec(i))).gt.0) then
      write(*,'(a)') trim(string_vec(i))
   end if
end do

end program rcsfmr


subroutine stringprocess(string1,string2)
implicit none
integer :: i,j,k,l,nstring,nunwanted,nints
character(len=500) :: string1,string2
character(len=1) :: unwanted(11),orb(9)

unwanted(1) = 'S'
unwanted(2) = 'P'
unwanted(3) = 'D'
unwanted(4) = 'F'
unwanted(5) = 'G'
unwanted(6) = 'H'
unwanted(7) = 'I'
unwanted(8) = 'K'
unwanted(9) = 'L'
unwanted(10) = '_'
unwanted(11) = '.'

orb(1) = 's'
orb(2) = 'p'
orb(3) = 'd'
orb(4) = 'f'
orb(5) = 'g'
orb(6) = 'h'
orb(7) = 'i'
orb(8) = 'k'
orb(9) = 'l'

string2 = ' '

string2(1:1) = string1(1:1)

! Remove unwanted strings

nstring = len(trim(string1))
k = 1
do i = 2,nstring
   nunwanted = 0
   do j = 1,11
      if (string1(i:i).eq.unwanted(j)) then
         nunwanted = 1
      end if
   end do
   if (nunwanted.eq.0) then
      k = k + 1
      string2(k:k) = string1(i:i)
   end if
end do

! Now remove digits unless digit in parenthesis or if
! first or second character to the right is 's','p',....

string1 = ' '
string1(1:1) = string2(1:1)
string1(2:2) = string2(2:2)

nstring = len(trim(string2))
k = 2
do i = 3,nstring
   nunwanted = 0
   select case (string2(i:i))
      case ('s')
         nunwanted = 1
      case ('p')
         nunwanted = 1
      case ('d')
         nunwanted = 1
      case ('f')
         nunwanted = 1
      case ('g')
         nunwanted = 1
      case ('h')
         nunwanted = 1
      case ('i')
         nunwanted = 1
      case ('k')
         nunwanted = 1
      case ('l')
         nunwanted = 1
      case ('(')
         nunwanted = 1
      case (')')
         nunwanted = 1
   end select
   if (nunwanted.eq.1) then
      k = k + 1
      string1(k:k) = string2(i:i)
   elseif (nunwanted.eq.0) then
      do j = 1,9
         if (string2(i+1:i+1).eq.orb(j)) nunwanted = 1
      end do
      if (string2(i-1:i-1).eq.'(') nunwanted = 1
      if (string2(i+1:i+1).eq.')') nunwanted = 1
      if (nunwanted.eq.1) then
         k = k + 1
         string1(k:k) = string2(i:i)
      end if
   end if

end do

!write(*,'(a)') trim(string1)

! Now add parenthesis and stars

string2 = ' '
string2(1:1) = string1(1:1)
string2(2:2) = string1(2:2)

nstring = len(trim(string1))
k = 2
do i = 3,nstring
   if (string1(i:i).eq.')') then
      k = k + 1
      string2(k:k) = ','
      k = k + 1
      string2(k:k) = '*'
      k = k + 1
      string2(k:k) = string1(i:i)
   else
      l = 0
      do j = 1,9
         if (string1(i:i).eq.orb(j)) then
            if (string1(i+1:i+1).ne.'(') then
               l = 1
               k = k + 1
               string2(k:k) = string1(i:i)
               k = k + 1
               string2(k:k) = '('
               k = k + 1
               string2(k:k) = '1'
               k = k + 1
               string2(k:k) = ','
               k = k + 1
               string2(k:k) = '*'
               k = k + 1
               string2(k:k) = ')'
            end if
         end if
      end do
      if (l.eq.0) then
         k = k + 1
         string2(k:k) = string1(i:i)
      end if
   end if
end do

!write(*,'(a)') trim(string2)



end
