program rseqtrans

! Per Jönsson, Malmö University, June 2015   

implicit none
logical :: ex
integer :: i,j,k,l,m,nstart,nstop,nstates,nplot,j2l,j2u
integer :: Zmin,Zmax,nions,Z(150),numplotu(100),numplotl(100)
integer :: posu,posl,nfound,readerr,nform,nfit,ntype
double precision :: etot,elev,esplit,transplot(100),A,gf,S
character(len=1) :: paru,parl,pplotu(100),pplotl(100),ans
character(len=2) :: mp,f1,f2,gauge
character(len=4) :: jvalu,jvall,jplotu(100),jplotl(100)
character(len=5) :: string
character(len=10) :: filename
character(len=34) :: mpstring,string2
character(len=100) :: label

write(*,*) 'RSEQTRANS'
write(*,*) 'This program reads output from rtransition for several'
write(*,*) 'ions and produces a Matlab/Octave file that plots'
write(*,*) 'A, gf, or S as a function of Z'  
write(*,*) 'Input files: transZ1, transZ2, .., transZn'
write(*,*) 'Output file: seqtransplot.m'
write(*,*) 

!--- Define Z range --------------------------

write(*,*) 'Give the first Z and last Z of the sequence'
read(*,*) Zmin,Zmax

!--- Give parity, J and number for state -----

write(*,*) 'Give multipolarity of transition: E1, M1, E2, M2'
read(*,*) mp

if (mp.eq.'E1') then
   mpstring =  ' Electric 2**( 1)-pole transitions'
elseif (mp.eq.'E2') then
   mpstring =  ' Electric 2**( 2)-pole transitions'
elseif (mp.eq.'M1') then
   mpstring =  ' Magnetic 2**( 1)-pole transitions'
elseif (mp.eq.'M2') then
   mpstring =  ' Magnetic 2**( 2)-pole transitions'
end if

write(*,*) 'How many transitions do you want to plot?'
read(*,*) nplot 
do i = 1,nplot
   write(*,*) 'Give number within symmetry,2*J and parity (+/-)'
   write(*,*) 'for upper and lower state'
   read(*,*) numplotu(i),j2u,pplotu(i),numplotl(i),j2l,pplotl(i)

   if (mod(j2u,2).eq.0) then 
      if (j2u.le.18) then
         write(jplotu(i),'(a3,i1)') '   ',j2u/2
      else
         write(jplotu(i),'(a2,i2)') '  ',j2u/2
      end if
   else
      if (j2u.le.9) then
         write(jplotu(i),'(a1,i1,a2)') ' ',j2u,'/2'
      else
         write(jplotu(i),'(i2,a2)') j2u,'/2'
      end if
   end if   
   if (mod(j2l,2).eq.0) then 
      if (j2l.le.18) then
         write(jplotl(i),'(a3,i1)') '   ',j2l/2
      else
         write(jplotl(i),'(a2,i2)') '  ',j2l/2
      end if
   else
      if (j2l.le.9) then
         write(jplotl(i),'(a1,i1,a2)') ' ',j2l,'/2'
      else
         write(jplotl(i),'(i2,a2)') j2l,'/2'
      end if
   end if   

!   write(*,*) pplotu(i),jplotu(i),numplotu(i)
!   write(*,*) pplotl(i),jplotl(i),numplotl(i)
end do

!--- What to plot and print -----------

write(*,*) 'Plot A (1), gf (2) or S (3) ?'
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

open(unit=12,file='seqtransplot.m',form='formatted')

!--- Some prepatory writing to the M-file ---

write(12,*) 'A = ['

j = 0
do i = Zmin,Zmax
   write(filename,'(a5,i0)') 'trans',i
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
!    Search for multipolarity string

   do
      read(100+k,'(a34)',iostat=readerr) string2
      if (readerr.ne.0) then
         write(*,*) 'The specified multipolarity not found'
         stop
      end if
      if (string2.eq.mpstring) then
         exit
      end if
   end do

   read(100+k,*)
   read(100+k,*)
   read(100+k,*)
   read(100+k,*)

!--- Now we are in the right position for reading --

   nfound = 0
   do
      read(100+k,300,iostat=readerr) f1,posu,jvalu,paru,f2,posl, &
!      read(100+k,300) f1,posu,jvalu,paru,f2,posl, &
                     jvall,parl,elev,gauge,A,gf,S 
      if (readerr.ne.0) then
         exit
      end if 
      if ((mp.eq.'E1').or.(mp.eq.'E2')) then
         read(100+k,301) gauge,A,gf,S            
      end if
!      write(*,300) f1,posu,jvalu,paru,f2,posl, &
!                     jvall,parl,elev,gauge,A,gf,S 
!      write(*,301) gauge,A,gf,S         

      do j = 1,nplot
         if ((pplotu(j).eq.paru).and.(jplotu(j).eq.jvalu).and.(numplotu(j).eq.posu).and. &
             (pplotl(j).eq.parl).and.(jplotl(j).eq.jvall).and.(numplotl(j).eq.posl)) then 
            if (nform.eq.1) then
               transplot(j) = A
            elseif (nform.eq.2) then
               transplot(j) = gf
            else 
               transplot(j) = S
            end if
            nfound = nfound + 1
         end if
      end do
   end do

98 continue   

   if (nfound.ne.nplot) then
      write(*,*) 'Specified states not found in all lists'
      stop
   end if
   write(12,*) Z(k),(transplot(j),j = 1,nplot)
end do

!--- Now finish the writing to the M-file ----

write(12,*) '];'
write(12,*) 'clf, hold on'
write(12,'(a,i3,a,i3,a)') 'zip = linspace(',Zmin,',',Zmax,');'
write(12,'(a)') "title('transition parameters as functions of Z')"
write(12,'(a)') "xlabel('Z')"
if (nform.eq.1) then
   write(12,'(a)') "ylabel('A (s-1)')"
elseif (nform.eq.2) then
   write(12,'(a)') "ylabel('gf')"
else
   write(12,'(a)') "ylabel('S')"
end if

write(12,*)
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
         write(12,'(a)') 'aiplsq = a(1)./zip.^2 + a(2)./zip + a(3) + a(4)*zip + a(5)*zip.^2 + a(6)*zip.^3;'
         write(12,'(a)') 'aiplsq = s*aiplsq + m;'
         write(12,'(a)') "plot(zip,aiplsq,'r')" 
         write(12,*)
      else
         write(12,'(a)') 'z = A(:,1);'
         write(12,'(a)') 'AD = [z.^0 z.^1 z.^2 z.^3];'
         write(12,'(a,i2,a)') 'y = A(:,',i+1,');'
         write(12,'(a)') 'm = mean(y); s = std(y);'
         write(12,'(a)') 'a = AD\(y-m)/s'         
         write(12,'(a)') 'aiplsq =  a(1) + a(2)*zip + a(3)*zip.^2 + a(4)*zip.^3;'
         write(12,'(a)') 'aiplsq = s*aiplsq + m;'
         write(12,'(a)') "plot(zip,aiplsq,'r')" 
         write(12,*)
      end if
   else
      write(12,'(a,i2,a)') "Aip = interp1(A(:,1),A(:,",i+1,"),zip,'spline');"
      write(12,*) 'plot(zip,Aip)'
      write(12,*)
   end if
end do

300 FORMAT(1X,A2,I3,1X,A4,1X,A1,2X,A2,I3,1X,A4,1X,A1,2X,0P,F13.2,A2,1P,3D13.5)
301 FORMAT(42X,A2,1P,3D13.5)

end program rseqtrans
   
