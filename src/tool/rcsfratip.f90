program rcsfratip

! Program to convert output format from grasp format to the 
! format used by RATIP

! Written by Jorgen Ekman and Per Jonsson, December 2013 

character(len=1000):: line1,line2,line3,line4,blankline
integer :: i,m,ncorr

do i = 1,1000
   blankline(i:i) = ' '
end do

open(unit=40,file='clist.out',status='old')
open(unit=41,file='clist_ratip.out',status='unknown')

do i = 1,5
  read(40,'(a)') line1
  write(41,'(a)') trim(line1)
end do

do
  read(40,'(a)',end=99) line1
  write(41,'(a)') trim(line1)
  read(40,'(a)') line2
  write(41,'(a)') trim(line2)

! Initialize new string

  line4 = blankline

! Analyze string 3

  read(40,'(a)') line3
  m = len_trim(line3)

  do
    do i = m,m-5,-1
      if (line3(i:i).eq.' ') then
        ncorr = 9*nint(real(i)/real(9))         ! Correct position of blank 
        line4(ncorr+1:ncorr+m-i) = line3(i+1:m) ! Move to correct position
        exit
       end if
    end do
    m = len_trim(line3(1:m-5))
    if (m.eq.0) exit
  end do
  write(41,'(a)') trim(line4)
end do

99 continue

end program rcsfratip
