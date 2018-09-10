program rseqenergy

! Per Jönsson, Malmö University, June 2015   

implicit none
logical :: ex
integer :: i,j,k,l,m,nstart,nstop,nstates,nplot,j2,nfit,ntype
integer :: Zmin,Zmax,nions,Z(150),numplot(100)
integer :: no,pos,nfound
double precision :: etot,elev,esplit,elevplot(100)
character(len=1) :: par,pplot(100),ans
character(len=4) :: jval,jplot(100)
character(len=5) :: string
character(len=10) :: filename
character(len=100) :: label

!

write(*,*) 'RSEQENERGY'
write(*,*) 'This program reads output from rlevels for several'
write(*,*) 'ions and produces a Matlab/Octave file that plots'
write(*,*) 'energy as a function of Z'  
write(*,*) 'Input files: energyZ1, energyZ2, .., energyZn'
write(*,*) 'Output file: seqenergyplot.m'
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
         write(jplot(i),'(a2,i1,a1)') '   ',j2/2,' '
      else
         write(jplot(i),'(a1,i2,a1)') '  ',j2/2,' '
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

open(unit=12,file='seqenergyplot.m',form='formatted')

!--- Some prepatory writing to the M-file ---

write(12,*) 'A = ['

j = 0
do i = Zmin,Zmax
   write(filename,'(a6,i0)') 'energy',i
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

!--- Find out where to start and stop reading ---

   i = 0
   j = 0
   do 
      read(100+k,'(a)',end=9) string
      j = j + 1
      if (string.eq.'-----') then
         i = i + 1
         if (i.eq.2) then 
            nstart = j
         end if
      end if
   end do

9  continue

   nstop = j

!   write(*,*) nstart,nstop
   nstates = nstop - nstart - 1

   rewind(unit=100+k)

!--- Read header information ---------------

   do i = 1,nstart
      read(100+k,*) 
   end do

!--- Read and xxxx

   nfound = 0
   do i = 1,nstates
      read(100+k,3) no,pos,jval,par,etot,elev,esplit,label
!      write(*,3) no,pos,jval,par,etot,elev,esplit,trim(label)
      do j = 1,nplot
         if ((pplot(j).eq.par).and.(jplot(j).eq.jval).and.(numplot(j).eq.pos)) then
            elevplot(j) = elev
            nfound = nfound + 1
         end if
      end do
   end do
   if (nfound.ne.nplot) then
      write(*,*) 'Specified states not found in all lists'
      stop
   end if
   write(12,*) Z(k),(elevplot(j),j = 1,nplot)
end do

!--- Now finish the writing to the M-file ----

write(12,*) '];'
write(12,*) 'clf, hold on'
write(12,'(a,i3,a,i3,a)') 'zip = linspace(',Zmin,',',Zmax,');'
do i = 1,nplot
   write(12,'(a,i2,a)') "plot(A(:,1),A(:,",i+1,"),'+')"
   write(12,'(a)') "title('E in cm-1 as a function of Z')"
   write(12,'(a)') "xlabel('Z')"
   write(12,'(a)') "ylabel('E (cm-1)')"
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

3  format(2i3,1x,a4,1x,a1,2x,f14.7,f12.2,f12.2,2x,a)   

end program rseqenergy
