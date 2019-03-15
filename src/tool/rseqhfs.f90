program rseqhfs

! Per Jönsson, Malmö University, June 2015

implicit none
logical :: ex
integer :: i,j,k,l,m,nstart,nstop,nstates,nplot,j2,nfit,ntype
integer :: Zmin,Zmax,nions,Z(150),numplot(100)
integer :: pos,nfound,nform,readerr
double precision :: A,B,gJ,hfsplot(100)
character(len=1) :: par,pplot(100),ans
character(len=4) :: jval,jplot(100)
character(len=5) :: string
character(len=10) :: filename
character(len=100) :: label,dummy

write(*,*) 'RSEQHFS'
write(*,*) 'This program reads output from rhfs for several'
write(*,*) 'ions and produces a Matlab/Octave file that '
write(*,*) 'plots hfs parameters as functions of Z'
write(*,*) 'Input files: hfsZ1, hfsZ2, .., hfsZn or'
write(*,*) 'Output file: seqhfsplot.m'
write(*,*)

!--- Define Z range --------------------------

write(*,*) 'Give the first Z and last Z of the sequence'
read(*,*) Zmin,Zmax

!--- Give parity, J and number for state -----

write(*,*) 'How many states do you want to plot?'
read(*,*) nplot
do i = 1,nplot
   write(*,*) 'Give number within symmetry,2*J and parity (+/-)'
   read(*,*) numplot(i),j2,pplot(i)

   if (mod(j2,2).eq.0) then
      if (j2.le.18) then
         write(jplot(i),'(a3,i1)') '   ',j2/2
      else
         write(jplot(i),'(a2,i2)') '  ',j2/2
      end if
   else
      if (j2.le.9) then
         write(jplot(i),'(a1,i1,a2)') ' ',j2,'/2'
      else
         write(jplot(i),'(i2,a2)') j2,'/2'
      end if
   end if
!   write(*,*) pplot(i),jplot(i),numplot(i)
end do

write(*,*) 'Plot A (1), B (2) or gJ (3) ?'
read(*,*) nform

write(*,*) 'Least-squares fit (y/n) ?'
read(*,*) ans
nfit = 0
if ((ans.eq.'y').or.(ans.eq.'Y')) then
   nfit = 1
   write(*,*) 'Type of fitting: a1 Z^-2 + a2 Z^-1 + ...+ a6 Z^3  (1)'
   write(*,*) '                 a1 + a2 Z + a3 Z^2 + a4 Z^3      (2)'
   read(*,*) ntype
end if

!--- Open existing files and define Z vector --

open(unit=12,file='seqhfsplot.m',form='formatted')

!--- Some prepatory writing to the M-file ---

write(12,*) 'A = ['

j = 0
do i = Zmin,Zmax
   write(filename,'(a3,i0)') 'hfs',i
   inquire(file=trim(filename),exist=ex)
   if (ex.eqv..true.) then
      j = j + 1
      Z(j) = i
      open(unit=100+j,file=trim(filename),form='formatted')
   end if
end do
nions = j

!write(*,*) 'nions',nions
!do i = 1,nions
!   write(*,*) Z(i)
!end do

!--- Loop over ions --------------------

do k = 1,nions

!--- Read header information ---------------

   do i = 1,9
      read(100+k,'(a)')
   end do

!--- Read and xxxx

   nfound = 0
   do
      read(100+k,3,iostat=readerr) pos,jval,par,A,B,gJ
       if (readerr.ne.0) then
         exit
      end if
!      write(*,3) pos,jval,par,A,B,gJ
      do j = 1,nplot
         if ((pplot(j).eq.par).and.(jplot(j).eq.jval).and.(numplot(j).eq.pos)) then
            if (nform.eq.1) then
               hfsplot(j) = A
            elseif (nform.eq.2) then
               hfsplot(j) = B
            else
               hfsplot(j) = gJ
            end if
            nfound = nfound + 1
         end if
      end do
   end do
   if (nfound.ne.nplot) then
      write(*,*) 'Specified states not found in all lists'
      stop
   end if
   write(12,*) Z(k),(hfsplot(j),j = 1,nplot)
end do

!--- Now finish the writing to the M-file ----

write(12,*) '];'
write(12,*) 'clf, hold on'
write(12,'(a,i3,a,i3,a)') 'zip = linspace(',Zmin,',',Zmax,');'
write(12,'(a)') "xlabel('Z')"
if (nform.eq.1) then
   write(12,'(a)') "ylabel('A_J (MHz)')"
elseif (nform.eq.2) then
   write(12,'(a)') "ylabel('B_J (MHz)')"
else
   write(12,'(a)') "ylabel('g_J')"
end if

do i = 1,nplot
   write(12,'(a,i2,a)') "plot(A(:,1),A(:,",i+1,"),'+')"
   write(12,*)
   if (nfit.eq.1) then
      if (ntype.eq.1) then
         write(12,'(a)') 'z = A(:,1);'
         write(12,'(a)') 'AD = [z.^(-2) z.^(-1) z.^0 z.^1 z.^2 z.^3];'
         write(12,'(a,i2,a)') 'y = A(:,',i+1,');'
         write(12,'(a)') 'm = mean(y); s = std(y);'
         write(12,'(a)') 'a = AD\(y-m)/s'
         write(12,'(a)') 'eiplsq = a(1)./zip.^2 + a(2)./zip + a(3) + a(4)*zip + a(5)*zip.^2 + a(6)*zip.^3;'
         write(12,'(a)') 'eiplsq = s*eiplsq + m;'
         write(12,'(a)') "plot(zip,eiplsq,'r')"
         write(12,*)
      else
         write(12,'(a)') 'z = A(:,1);'
         write(12,'(a)') 'AD = [z.^0 z.^1 z.^2 z.^3];'
         write(12,'(a,i2,a)') 'y = A(:,',i+1,');'
         write(12,'(a)') 'm = mean(y); s = std(y);'
         write(12,'(a)') 'a = AD\(y-m)/s'
         write(12,'(a)') 'eiplsq =  a(1) + a(2)*zip + a(3)*zip.^2 + a(4)*zip.^3;'
         write(12,'(a)') 'eiplsq = s*eiplsq + m;'
         write(12,'(a)') "plot(zip,eiplsq,'r')"
         write(12,*)
      end if
   else
      write(12,'(a,i2,a)') "eip = interp1(A(:,1),A(:,",i+1,"),zip,'spline');"
      write(12,*) 'plot(zip,eip)'
      write(12,*)
   end if
end do

3  format(i4,5x,a4,1x,a1,2x,1P,3D20.10)

end program rseqhfs
