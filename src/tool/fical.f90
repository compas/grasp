program fical
  implicit none
  integer :: i,j,ii,jj,l,lsjpresent,lmax,lmax1,lmax2
  integer :: nvec1,nvec2,nveclsj1,nveclsj2,level1(500),level2(500),levellsj1(500),levellsj2(500)
  integer :: relflag,dspin,efacflag

  double precision :: energy1(500),nms1(500),sms1(500),ms1(500),fs1(6,500),fshift1(500),fshift11(500),fspin1(500)
  double precision :: trms_1,t_1,w_1,b20_1,b40_1, radmom1(4)
  double precision :: trms_2,t_2,w_2,b20_2,b40_2, radmom2(4)

  double precision :: energy2(500),nms2(500),sms2(500),ms2(500),fs2(6,500),fshift2(500),fshift22(500),fspin2(500)
  double precision :: tenergy, nmslaw, nmsconst, m1, m2, dr2, masshift,fieldshift,fieldshift22,fieldshift1(500),fieldshift2(500)
  double precision :: isotopeshift
  double precision :: f0, f2, f4, f6, f0res, f2res, linfac
  double precision :: fieldshiftres, fieldshift22res, linfac2

  double precision :: drms, au2Kays
  double precision :: z

  character(len=100) :: state1, state2, file1, file2, file3, file4, file5, file6
  character(len=11) :: spinpar1(500),spinpar2(500)
  character(len=4) :: typeflag,spin1(500),spin2(500),spinlsj1(500),spinlsj2(500)
  character(len=3) :: unit
  character(len=1) :: def,ci,rel,efac,corr,soph,par1(500),par2(500),parlsj1(500),parlsj2(500)
  character(len=200) :: csflabel1(500),csflabel2(500),csflabel1scal,csflabel2scal
  character(len=200) :: headermpGHz,headermpMHz,headermpmeV
  character(len=200) :: headermpGHzlsj,headermpMHzlsj,headermpmeVlsj
  character(len=200) :: headerfpGHz,headerfpMHz,headerfpmeV
  character(len=200) :: headerfpGHzlsj,headerfpMHzlsj,headerfpmeVlsj
  character(len=200) :: headerfp2GHz,headerfp2MHz,headerfp2meV
  character(len=200) :: headerfp2GHzlsj,headerfp2MHzlsj,headerfp2mevlsj
  character(len=200) :: headereGHz,headereMHz,headeremeV
  character(len=200) :: headereGHzlsj,headereMHzlsj,headeremeVlsj
  character(len=50) :: format1a,format1b,format1b2,format1c,format1c2,format2
  character(len=50) :: format1alsj,format1blsj,format1b2lsj,format1clsj,format1c2lsj,format2lsj
  character(len=200) :: uplevel,lowlevel

  logical :: exist1,exist2
  parameter (au2Kays = 219474.63136320d0)   ! CODATA 2018

  headermpGHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (GHz u)" &
       ,6X,"NMS (GHz u)",6X,"SMS (GHz u)"8X,"MS (GHz u)")'
  headermpGHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (GHz u)" &
       ,6X,"NMS (GHz u)",6X,"SMS (GHz u)"8X,"MS (GHz u)")'
  headermpMHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (MHz u)" &
       ,6X,"NMS (MHz u)",6X,"SMS (MHz u)"8X,"MS (MHz u)")'
  headermpMHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (MHz u)" &
       ,6X,"NMS (MHz u)",6X,"SMS (MHz u)"8X,"MS (MHz u)")'
  headermpmeV = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",3X,"NMS-S (meV u)" &
       ,6X,"NMS (meV u)",6X,"SMS (meV u)"8X,"MS (meV u)")'
  headermpmeVlsj = '(X,a,5X,a,5X,"Energy (cm-1)",3X,"NMS-S (meV u)" &
       ,6X,"NMS (meV u)",6X,"SMS (meV u)"8X,"MS (meV u)")'

  headerfpGHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,3X,"F0 (GHz fm-2)",5X,"F2 (GHz fm-4)",5X,"F4 (GHz fm-6)",5X,"F6 (GHz fm-8)",5X,"Lin. dep. (per mille)")'
  headerfpGHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,3X,"F0 (GHz fm-2)",5X,"F2 (GHz fm-4)",5X,"F4 (GHz fm-6)",5X,"F6 (GHz fm-8)",5X,"Lin. dep. (per mille)")'
  headerfpMHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,3X,"F0 (MHz fm-2)",5X,"F2 (MHz fm-4)",5X,"F4 (MHz fm-6)",5X,"F6 (MHz fm-8)")'
  headerfpMHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,3X,"F0 (MHz fm-2)",5X,"F2 (MHz fm-4)",5X,"F4 (MHz fm-6)",5X,"F6 (MHz fm-8)")'
  headerfpmeV = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,3X,"F0 (meV fm-2)",5X,"F2 (meV fm-4)",5X,"F4 (meV fm-6)",5X,"F6 (meV fm-8)")'
  headerfpmeVlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,3X,"F0 (meV fm-2)",5X,"F2 (meV fm-4)",5X,"F4 (meV fm-6)",5X,"F6 (meV fm-8)")'


  headerfp2GHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (GHz fm-2)",1X,"F0VED1 (GHz fm-4)",3X,"Lin. dep. (per mille)")'
  headerfp2GHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (GHz fm-2)",1X,"F0VED1 (GHz fm-4)",3X,"Lin. dep. (per mille)")'
  headerfp2MHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (MHz fm-2)",1X,"F0VED1 (MHz fm-4)")'
  headerfp2MHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (MHz fm-2)",1X,"F0VED1 (MHz fm-4)")'
  headerfp2meV = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (meV fm-2)",1X,"F0VED1 (meV fm-4)")'
  headerfp2meVlsj = '(X,a,5X,a,5X,"Energy (cm-1)" &
       ,1X,"F0VED0 (meV fm-2)",1X,"F0VED1 (meV fm-4)")'


  headereGHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (GHz)",9X,"FS (GHz)",10X,"IS (GHz)")'
  headereGHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (GHz)",9X,"FS (GHz)",10X,"IS (GHz)")'
  headereMHz = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (MHz)",9X,"FS (MHz)",10X,"IS (MHz)")'
  headereMHzlsj = '(X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (MHz)",9X,"FS (MHz)",10X,"IS (MHz)")'
  headeremeV = '(X,"Upper level",5X,"Lower level",5X,"Energy (cm-1)",6X,"MS (meV)",9X,"FS (meV)",10X,"IS (meV)")'
  headeremeVlsj = '(X,a,5X,a,5X,"Energy (cm-1)",6X,"MS (meV)",9X,"FS (meV)",10X,"IS (meV)")'

!  format1a = '(i4,a13,i3,a12,f12.2,1P,4d16.7)'
!  format1b = '(i4,a13,i3,a12,f12.2,1P,4d16.7)'
!  format1c = '(i4,a13,i3,a12,f12.2,1P,2d16.7)'
!  format2 = '(i4,a13,i3,a12,f12.2,1P,3d16.7)'

  format1a = '(i4,a13,i3,a12,f12.2,1P,4d18.7)'
  format1alsj = '(a,a10,a,a10,f12.2,1P,4d18.7)'
  format1b = '(i4,a13,i3,a12,f12.2,1P,4d18.7)'
  format1blsj = '(a,a10,a,a10,f12.2,1P,4d18.7)'
  format1b2 = '(i4,a13,i3,a12,f12.2,1P,4d18.7,0P,f18.7)'
  format1b2lsj = '(a,a10,a,a10,f12.2,1P,4d18.7,0P,f18.7)'

  format1c = '(i4,a13,i3,a12,f12.2,1P,2d18.7)'
  format1clsj = '(a,a10,a,a10,f12.2,1P,2d18.7)'
  format1c2 = '(i4,a13,i3,a12,f12.2,1P,2d18.7,0P,f18.7)'
  format1c2lsj = '(a,a10,a,a10,f12.2,1P,2d18.7,0P,f18.7)'

  format2 = '(i4,a13,i3,a12,f12.2,1P,3d18.7)'
  format2lsj = '(a,a10,a,a10,f12.2,1P,3d18.7)'

  write(*,*) 'WELCOME TO PROGRAM FICAL'
  write(*,*) 'Computes line frequency isotope shift parameters and/or energies using'
  write(*,*) 'output files from ris4'
  write(*,*)
  write(*,*) 'Input files: <state1>.(c)i, <state2>.(c)i'
  write(*,*) 'Output file: <state1>.<state2>.(c)fi'

  ! READ IN USER DEFINED PARAMETERS
  write(*,*)
  write(*,*) 'Default settings (y/n)?: '
  read(*,*) def
  if(def.ne.'y'.and.def.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if

  write(*,*) 'Give name of state 1: '
  read(*,'(a100)') state1
  write(*,*) 'Give name of state 2: '
  read(*,'(a100)') state2

  write(*,*) 'Resulting isotope shifts from CI calculations (y/n)?: '
  read(*,*) ci
  if(ci.ne.'y'.and.ci.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if

  write(*,*) 'Have electronic factors been calculated (y/n)?: '
  read(*,*) efac
  if(efac.ne.'y'.and.efac.ne.'n') then
     write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
     call exit()
  end if

  if(def.eq.'y') then
     typeflag = 'both'
     unit  = 'GHz'
     rel = 'y'
     soph = 'n'
  else if(def.eq.'n') then
     write(*,*) 'Compute IS parameters (para), IS energies (ener) or both (both)?: '
     read(*,*) typeflag
     if(typeflag.ne.'para'.and.typeflag.ne.'ener'.and.typeflag.ne.'both') then
        write(*,*) 'Your input must be "para", "ener" or "both" (case sensitive). Try again!'
        call exit()
     end if

     write(*,*) 'Units (GHz, MHz or meV)?: '
     read(*,*) unit
     if(unit.ne.'GHz'.and.unit.ne.'MHz'.and.unit.ne.'meV') then
        write(*,*) 'Your input must be "GHz", "MHz", or "meV" (case sensitive). Try again!'
        call exit()
     end if

     write(*,*) 'Use relativistically corrected mass shift parameters (y/n)?: '
     read(*,*) rel
     if(rel.ne.'y'.and.rel.ne.'n') then
        write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
        call exit()
     end if

     if(typeflag.eq.'ener'.or.typeflag.eq.'both') then
        write(*,*) 'Use sophisticated model for radial moments (y/n)?: '
        read(*,*) soph
        if(soph.ne.'y'.and.soph.ne.'n') then
           write(*,*) 'Your input must be "y" or "n" (case sensitive). Try again!'
           call exit()
        end if
     else if(typeflag.eq.'para') then
        soph = 'n'
     end if
  end if

  if ((typeflag.eq.'ener').or.(typeflag.eq.'both')) then
     if(soph.eq.'y') then
        write(*,*) 'Data for isotope 1'
        write(*,*) 'Enter mass(amu),rms radius, t, omega, b20, b40: '
        read(*,*) m1,trms_1,t_1,w_1,b20_1,b40_1
        write(*,*) 'Data for isotope 2'
        write(*,*) 'Enter mass(amu),rms radius, t, omega, b20, b40: '
        read(*,*) m2,trms_2,t_2,w_2,b20_2,b40_2
     else if(soph.eq.'n') then
        write(*,*) 'Mass of isotope 1 (in amu)?: '
        read(*,*) m1
        write(*,*) 'Mass of isotope 2 (in amu)?: '
        read(*,*) m2
        write(*,*) 'd<r^2> (isotope 1 - isotope 2) (in fm^2)?: '
        read(*,*) dr2
     end if
  end if

  ! NAME INPUT AND OUTPUT FILES

  file4 = 'isodata'
  l = index(state1,' ')
  file5 = state1(1:l-1)//'.lsj.lbl'
  l = index(state2,' ')
  file6 = state2(1:l-1)//'.lsj.lbl'
  if(ci.eq.'y') then
     l = index(state1,' ')
     file1 = state1(1:l-1)//'.ci'
     l = index(state2,' ')
     file2 = state2(1:l-1)//'.ci'
     file3 = trim(state1)//'.'//trim(state2)//'.cfi'
     l = index(state1,' ')
  else
     l = index(state1,' ')
     file1 = state1(1:l-1)//'.i'
     l = index(state2,' ')
     file2 = state2(1:l-1)//'.i'
     file3 = trim(state1)//'.'//trim(state2)//'.fi'
     l = index(state1,' ')
  end if

  ! DEFINE CERTAIN FLAGS

  if(rel.eq.'y') then
     relflag = 1
  else
     relflag = 0
  end if

  if(efac.eq.'y') then
     efacflag = 1
  else
     efacflag = 0
  end if

  ! OPEN AND READ IN DATA FROM INPUT FILES

  lsjpresent = 0
  inquire(file=file5, exist=exist1)
  inquire(file=file6, exist=exist2)
  if(exist1.and.exist2) then
     lsjpresent = 1
     call openlsj(file5,nveclsj1,levellsj1,spinlsj1,parlsj1,csflabel1,lmax1)
     call openlsj(file6,nveclsj2,levellsj2,spinlsj2,parlsj2,csflabel2,lmax2)
     lmax = max(lmax1,lmax2)
     uplevel = "Upper level                                                   "
     lowlevel = "Lower level                                                   "

  end if

  call openfile(file1,unit,efacflag,relflag,nvec1,level1,spinpar1,fspin1,par1,energy1,nms1,sms1,ms1,fs1)

  call openfile(file2,unit,efacflag,relflag,nvec2,level2,spinpar2,fspin2,par2,energy2,nms2,sms2,ms2,fs2)

  call openisodata(file4,file3,z)

  ! WRITE INPUT ISOTOPE DATA TO OUPTPUT FILE
  if(typeflag.ne.'para') then
     write(7,'(a)') ' INPUT ISOTOPE DATA'
     if(soph.eq.'n') then
        write(7,'(a22,f10.6,a)') ' Mass of isotope 1: ',m1,' amu'
        write(7,'(a22,f10.6,a)') ' Mass of isotope 2: ',m2,' amu'
        write(7,'(a22,f10.6,a)') ' d<r^2>: ',dr2,' fm^2'
     else
        write(7,*) '               Isotope 1    Isotope 2'
        write(7,'(a16,f10.6,f13.6)') ' Mass [amu]: ',m1,m2
        write(7,'(a16,f10.6,f13.6)') ' r_rms [fm]: ',trms_1,trms_2
        write(7,'(a16,f10.6,f13.6)') ' t [fm]: ',t_1,t_2
        write(7,'(a16,f10.6,f13.6)') ' omega: ',w_1,w_2
        write(7,'(a16,f10.6,f13.6)') ' b_20: ',b20_1,b20_2
        write(7,'(a16,f10.6,f13.6)') ' b_40: ',b40_1,b40_2
     end if
  end if
  write(7,*)

  if(efacflag.eq.1) then
     fshift1(:) = fs1(5,:)
     fshift2(:) = fs2(5,:)
     fshift11(:) = fs1(6,:)
     fshift22(:) = fs2(6,:)
  else
     if(unit.eq.'GHz') then
        fs1(1,:) = fs1(1,:)*z*4.9210908105D-3
        fs2(1,:) = fs2(1,:)*z*4.9210908105D-3
     else if(unit.eq.'MHz') then
        fs1(1,:) = fs1(1,:)*z*4.9210908105d0
        fs2(1,:) = fs2(1,:)*z*4.9210908105d0
     else if(unit.eq.'meV') then
        fs1(1,:) = fs1(1,:)*z*4.9210908105D-3/241.7989348d0
        fs2(1,:) = fs2(1,:)*z*4.9210908105D-3/241.7989348d0
     end if
     fshift1(:) = fs1(1,:)
     fshift2(:) = fs2(1,:)
     fshift11(:) = fs1(1,:)
     fshift22(:) = fs2(1,:)
  end if

  if(soph.eq.'y') then
     call fermidist(trms_1,t_1,w_1,b20_1,b40_1,0, radmom1)
     call fermidist(trms_2,t_2,w_2,b20_2,b40_2,0, radmom2)

     drms =abs(radmom1(1)-radmom2(1))

     fieldshift1(:) = 0.0d0
     do i=1,nvec1
        do j=1,4
           fieldshift1(i) = fieldshift1(i) + fs1(j,i)*(radmom1(j)-radmom2(j))
        end do
     end do
     fieldshift2(:) = 0.0d0
     do i=1,nvec2
        do j=1,4
           fieldshift2(i) = fieldshift2(i) + fs2(j,i)*(radmom1(j)-radmom2(j))
        end do
     end do

     ! WRITE NUCLEAR RADIAL MOMENTS TO OUPTPUT FILE
     if(typeflag.ne.'para') then
        write(7,'(a)') ' NUCLEAR RADIAL MOMENTS'
        write(7,*) '               Isotope 1    Isotope 2     Isotope 1 - Isotope 2'
        write(7,'(a16,1P,D11.5,D13.5,D15.5)') ' <r^2> [fm^2]: ',radmom1(1),radmom2(1),radmom1(1)-radmom2(1)
        write(7,'(a16,1P,D11.5,D13.5,D15.5)') ' <r^4> [fm^4]: ',radmom1(2),radmom2(2),radmom1(2)-radmom2(2)
        write(7,'(a16,1P,D11.5,D13.5,D15.5)') ' <r^6> [fm^6]: ',radmom1(3),radmom2(3),radmom1(3)-radmom2(3)
        write(7,'(a16,1P,D11.5,D13.5,D15.5)') ' <r^8> [fm^8]: ',radmom1(4),radmom2(4),radmom1(4)-radmom2(4)
        write(7,*)
     end if

  end if


  ! COMPUTE AND WRITE LINE ISOTOPE SHIFT PARAMETERS TO FILE

  if((typeflag.eq.'para').or.(typeflag.eq.'both')) then
     write(7,*) 'LINE MASS SHIFT PARAMETERS'
     write(7,*)
     if(unit.eq.'GHz') then
        nmsconst = 29.9792458d0/1822.888485d0
        if(lsjpresent.eq.0) then
           write(7,headermpGHz)
        else
           write(7,headermpGHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'MHz') then
        nmsconst = 29979.2458d0/1822.888485d0
        if(lsjpresent.eq.0) then
           write(7,headermpMHz)
        else
           write(7,headermpMHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'meV') then
        nmsconst = 0.12398419739d0/1822.888485d0
        if(lsjpresent.eq.0) then
           write(7,headermpmeV)
        else
           write(7,headermpmeVlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     end if
     write(7,*)
     do i=1,nvec1
        do ii=1,nveclsj1
           if((level1(i).eq.levellsj1(ii)).and.(spinpar1(i)(6:9).eq.spinlsj1(ii)).and.(spinpar1(i)(11:11).eq.parlsj1(ii))) then
              csflabel1scal = csflabel1(ii)
           end if
        end do
        do j=1,nvec2
           do jj=1,nveclsj2
              if((level2(j).eq.levellsj2(jj)).and.(spinpar2(j)(6:9).eq.spinlsj2(jj)).and.(spinpar2(j)(11:11).eq.parlsj2(jj))) then
                 csflabel2scal = csflabel2(jj)
              end if
           end do
           if(energy1(i).gt.energy2(j)) then
              tenergy = (energy1(i)-energy2(j))*Au2kays
              nmslaw = -1.d0*tenergy*nmsconst
              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1a) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                 else
                    write(7,format1alsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1a) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                 else
                    write(7,format1alsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,nmslaw,nms1(i)-nms2(j),sms1(i)-sms2(j),ms1(i)-ms2(j)
                 end if
              end if
           else
              tenergy = (energy2(j)-energy1(i))*au2Kays
              nmslaw = -1.d0*tenergy*nmsconst
              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1a) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                 else
                    write(7,format1alsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1a) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                 else
                    write(7,format1alsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,nmslaw,nms2(j)-nms1(i),sms2(j)-sms1(i),ms2(j)-ms1(i)
                 end if
              end if
           end if
        end do
     end do
!     if(efacflag.eq.1) then
     write(7,*)
     write(7,*) 'LINE FIELD SHIFT PARAMETERS'
     write(7,*)
     if(unit.eq.'GHz') then
        if(lsjpresent.eq.0) then
           write(7,headerfpGHz)
        else
           write(7,headerfpGHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'MHz') then
        if(lsjpresent.eq.0) then
           write(7,headerfpMHz)
        else
           write(7,headerfpMHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'meV') then
        if(lsjpresent.eq.0) then
           write(7,headerfpmeV)
        else
           write(7,headerfpmeVlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     end if
     write(7,*)
     do i=1,nvec1
        do ii=1,nveclsj1
           if((level1(i).eq.levellsj1(ii)).and.(spinpar1(i)(6:9).eq.spinlsj1(ii)).and.(spinpar1(i)(11:11).eq.parlsj1(ii))) then
              csflabel1scal = csflabel1(ii)
           end if
        end do
        do j=1,nvec2
           do jj=1,nveclsj2
              if((level2(j).eq.levellsj2(jj)).and.(spinpar2(j)(6:9).eq.spinlsj2(jj)).and.(spinpar2(j)(11:11).eq.parlsj2(jj))) then
                 csflabel2scal = csflabel2(jj)
              end if
           end do
           if(energy1(i).gt.energy2(j)) then
              tenergy = (energy1(i)-energy2(j))*au2Kays
              f0 = fs1(1,i)-fs2(1,j)
              f2 = fs1(2,i)-fs2(2,j)
              f4 = fs1(3,i)-fs2(3,j)
              f6 = fs1(4,i)-fs2(4,j)
              fieldshift = fshift1(i)-fshift2(j)
              fieldshift22 = fshift11(i)-fshift22(j)
              if(i.eq.1.and.j.eq.1) then
                 f0res = f0
                 f2res = f2
              end if
              linfac = 1000.d0*(f0*f2res/(f2*f0res)-1.d0)

              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1b) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,f0,f2,f4,f6
                 else
                    write(7,format1blsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,f0,f2,f4,f6
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1b2) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,f0,f2,f4,f6,linfac
                 else
                    write(7,format1b2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,f0,f2,f4,f6,linfac
                 end if
              end if
           else
              tenergy = (energy2(j)-energy1(i))*au2Kays
              f0 = fs2(1,j)-fs1(1,i)
              f2 = fs2(2,j)-fs1(2,i)
              f4 = fs2(3,j)-fs1(3,i)
              f6 = fs2(4,j)-fs1(4,i)
              fieldshift = fshift2(j)-fshift1(i)
              fieldshift22 = fshift22(j)-fshift11(i)
              if(i.eq.1.and.j.eq.1) then
                 f0res = f0
                 f2res = f2
              end if
              linfac = 1000.d0*(f0*f2res/(f2*f0res)-1.d0)

              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1b) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,f0,f2,f4,f6
                 else
                    write(7,format1blsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,f0,f2,f4,f6
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1b2) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,f0,f2,f4,f6,linfac
                 else
                    write(7,format1b2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,f0,f2,f4,f6,linfac
                 end if

              end if
           end if
        end do
     end do

     ! end if
  !end if

!     write(7,*)
!     write(7,*) 'LINE FIELD SHIFT PARAMETERS 2' ! JE
     write(7,*)
     if(unit.eq.'GHz') then
        if(lsjpresent.eq.0) then
           write(7,headerfp2GHz)
        else
           write(7,headerfp2GHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'MHz') then
        if(lsjpresent.eq.0) then
           write(7,headerfp2MHz)
        else
           write(7,headerfp2MHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'meV') then
        if(lsjpresent.eq.0) then
           write(7,headerfp2meV)
        else
           write(7,headerfp2meVlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     end if
     write(7,*)
     do i=1,nvec1
        do ii=1,nveclsj1
           if((level1(i).eq.levellsj1(ii)).and.(spinpar1(i)(6:9).eq.spinlsj1(ii)).and.(spinpar1(i)(11:11).eq.parlsj1(ii))) then
              csflabel1scal = csflabel1(ii)
           end if
        end do
        do j=1,nvec2
           do jj=1,nveclsj2
              if((level2(j).eq.levellsj2(jj)).and.(spinpar2(j)(6:9).eq.spinlsj2(jj)).and.(spinpar2(j)(11:11).eq.parlsj2(jj))) then
                 csflabel2scal = csflabel2(jj)
              end if
           end do
           if(energy1(i).gt.energy2(j)) then
              tenergy = (energy1(i)-energy2(j))*au2Kays
              f0 = fs1(1,i)-fs2(1,j)
              f2 = fs1(2,i)-fs2(2,j)
              f4 = fs1(3,i)-fs2(3,j)
              f6 = fs1(4,i)-fs2(4,j)
              fieldshift = fshift1(i)-fshift2(j)
              fieldshift22 = fshift11(i)-fshift22(j)
              if(i.eq.1.and.j.eq.1) then
                 fieldshiftres = fieldshift
                 fieldshift22res = fieldshift22
              end if
              linfac2 = 1000.d0*(fieldshift*fieldshift22res/(fieldshift22*fieldshiftres)-1.d0)

              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1c) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,fieldshift,fieldshift22
                 else
                    write(7,format1clsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,fieldshift,fieldshift22
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1c2) level1(i), adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))), &
                         tenergy,fieldshift,fieldshift22,linfac2
                 else
                    write(7,format1c2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,fieldshift,fieldshift22,linfac2
                 end if
              end if
           else
              tenergy = (energy2(j)-energy1(i))*au2Kays
              f0 = fs2(1,j)-fs1(1,i)
              f2 = fs2(2,j)-fs1(2,i)
              f4 = fs2(3,j)-fs1(3,i)
              f6 = fs2(4,j)-fs1(4,i)
              fieldshift = fshift2(j)-fshift1(i)
              fieldshift22 = fshift22(j)-fshift11(i)
              if(i.eq.1.and.j.eq.1) then
                 fieldshiftres = fieldshift
                 fieldshift22res = fieldshift22
              end if
              linfac2 = 1000.d0*(fieldshift*fieldshift22res/(fieldshift22*fieldshiftres)-1.d0)

              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format1c) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,fieldshift,fieldshift22
                 else
                    write(7,format1clsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,fieldshift,fieldshift22
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format1c2) level2(j), adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))), &
                         tenergy,fieldshift,fieldshift22,linfac2
                 else
                    write(7,format1c2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,fieldshift,fieldshift22,linfac2
                 end if
              end if
           end if
        end do
     end do

     ! end if
  end if

  ! COMPUTE AND WRITE LINE ISOTOPE SHIFT ENERGIES TO FILE

  write(7,*)
  if((typeflag.eq.'ener').or.(typeflag.eq.'both')) then
     write(7,*) 'LINE ISOTOPE SHIFT ENERGIES'
     write(7,*)
     if(unit.eq.'GHz') then
        if(lsjpresent.eq.0) then
           write(7,headereGHz)
        else
           write(7,headereGHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'MHz') then
        if(lsjpresent.eq.0) then
           write(7,headereMHz)
        else
           write(7,headereMHzlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     else if(unit.eq.'meV') then
        if(lsjpresent.eq.0) then
           write(7,headeremeV)
        else
           write(7,headeremeVlsj) uplevel(1:lmax+8),lowlevel(1:lmax+8)
        end if
     end if

     write(7,*)
     do i=1,nvec1
        do ii=1,nveclsj1
           if((level1(i).eq.levellsj1(ii)).and.(spinpar1(i)(6:9).eq.spinlsj1(ii)).and.(spinpar1(i)(11:11).eq.parlsj1(ii))) then
              csflabel1scal = csflabel1(ii)
           end if
        end do
        do j=1,nvec2
           do jj=1,nveclsj2
              if((level2(j).eq.levellsj2(jj)).and.(spinpar2(j)(6:9).eq.spinlsj2(jj)).and.(spinpar2(j)(11:11).eq.parlsj2(jj))) then
                 csflabel2scal = csflabel2(jj)
              end if
           end do
           if(energy1(i).gt.energy2(j)) then
              tenergy = (energy1(i)-energy2(j))*au2Kays
              if(m1.gt.m2) then
                 masshift = (ms1(i)-ms2(j))*(m2-m1)/(m1*m2)
                 if(soph.eq.'y') then
                    fieldshift = fieldshift1(i)-fieldshift2(j)
                 else if(soph.eq.'n') then
                    fieldshift = (fshift1(i)-fshift2(j))*dr2 + (fshift11(i)-fshift22(j))*dr2**2.d0
                 end if
              else
                 masshift = (ms1(i)-ms2(j))*(m1-m2)/(m1*m2)
                 if(soph.eq.'y') then
                    fieldshift = -1.d0*(fieldshift1(i)-fieldshift2(j))
                 else if(soph.eq.'n') then
                    fieldshift = -1.d0*(fshift2(j)-fshift1(i))*dr2 + (fshift22(j)-fshift11(i))*dr2**2.d0
                 end if
              end if
              isotopeshift = masshift + fieldshift
              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format2) level1(i),adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))),tenergy, &
                         masshift,fieldshift,isotopeshift
                 else
                    write(7,format2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,masshift,fieldshift,isotopeshift
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format2) level1(i),adjustl(trim(spinpar1(i))),level2(j),adjustl(trim(spinpar2(j))),tenergy, &
                         masshift,fieldshift,isotopeshift
                 else
                    write(7,format2lsj) csflabel1scal(1:lmax+3),adjustl(trim(spinpar1(i))),csflabel2scal(1:lmax+3), &
                         adjustl(trim(spinpar2(j))),tenergy,masshift,fieldshift,isotopeshift
                 end if
              end if
           else
              tenergy = (energy2(j)-energy1(i))*au2Kays
              if(m1.gt.m2) then
                 masshift = (ms2(j)-ms1(i))*(m2-m1)/(m1*m2)
                 if(soph.eq.'y') then
                    fieldshift = fieldshift2(j)-fieldshift1(i)
                 else if(soph.eq.'n') then
                    !write(*,*) 'm1 larger than m2'
                    fieldshift = (fshift2(j)-fshift1(i))*dr2 + (fshift22(j)-fshift11(i))*dr2**2.d0
                 end if
              else
                 masshift = (ms2(j)-ms1(i))*(m1-m2)/(m1*m2)
                 if(soph.eq.'y') then
                    fieldshift = -1.d0*(fieldshift2(j)-fieldshift1(i))
                 else if(soph.eq.'n') then
                    !write(*,*) 'm1 smaller than m2'
                    !write(*,*) fshift2(j)-fshift1(i),fshift22(j)-fshift11(i)
                    fieldshift = -1.d0*(fshift2(j)-fshift1(i))*dr2 + (fshift22(j)-fshift11(i))*dr2**2.d0
                 end if
              end if
              isotopeshift = masshift + fieldshift
              if(state1.eq.state2.and.j.gt.i) then
                 if(lsjpresent.eq.0) then
                    write(7,format2) level2(j),adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))),tenergy, &
                         masshift,fieldshift,isotopeshift
                 else
                    write(7,format2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,masshift,fieldshift,isotopeshift
                 end if
              else if(state1.ne.state2) then
                 if(lsjpresent.eq.0) then
                    write(7,format2) level2(j),adjustl(trim(spinpar2(j))),level1(i),adjustl(trim(spinpar1(i))),tenergy, &
                         masshift,fieldshift,isotopeshift
                 else
                    write(7,format2lsj) csflabel2scal(1:lmax+3),adjustl(trim(spinpar2(j))),csflabel1scal(1:lmax+3), &
                         adjustl(trim(spinpar1(i))),tenergy,masshift,fieldshift,isotopeshift
                 end if
              end if
           end if
        end do
     end do
  end if
  write(7,*)
  close(7)
  write(*,*)
  write(*,*) 'program FICAL finished ...'
  write(*,*)
  write(*,*) 'Isotope shift parameters/energies written to file ', file3
end program fical

subroutine openisodata(name,nameout,z)
  implicit none
  character(len=200) :: string
  character(len=100) :: name, nameout
  character(len=21) :: dummy21
  character(len=21) :: dummy6
  integer :: i
  double precision :: z,a,c,t,rms,bohr,pi,radmom3(4)

  bohr = 52917.721067d0
  pi = 3.14159265359d0

  open (7,file=nameout,status='unknown',form='formatted')
  open (8,file=name,status='unknown',form='formatted')

  read(8,'(a)') string
  read(8,'(d18.5)') z
  read(8,'(a)') string
  read(8,'(a)') string
  read(8,'(a)') string
  read(8,'(d18.5)') a
  read(8,'(a)') string
  read(8,'(d18.5)') c

  a = a*4.d0*log(3.d0)
  call fermidist(c,a,0.0d0,0.0d0,0.0d0,1, radmom3)

  write(7,*)
  write(7,'(a)') " REFERENCE ISOTOPE DATA FROM ISODATA"
  write(7,'(a,f10.6)') " Atomic number: ", z
  write(7,'(a)') " Fermi nucleus: "
  write(7,'(a16,f10.6,a)') "c: ", c, " fm"
  write(7,'(a16,f10.6,a)') "r_rms: ", sqrt(radmom3(1)), " fm"
  write(7,'(a16,f10.6,a)') "t: ", a, " fm"
  write(7,*)

  close(8)

end subroutine openisodata

subroutine openfile(name,unit,eflag,rflag,nvec,level,spinpar,fspin,par,energy,nms,sms,ms,fs)
  implicit none
  integer :: i,j,l
  integer :: nvec,level(500)
  integer :: eflag,rflag
  double precision :: energy(500),nms(500),nmscorr(500),nmsrel(500)
  double precision :: sms(500),smscorr(500),smsrel(500)
  double precision :: ms(500),mscorr(500),msrel(500)
  double precision :: fs(6,500),fspin(500)
  double precision :: meV2GHz
  character(len=200) :: string
  character(len=100) :: name
  character(len=24) :: dummy24
  character(len=23) :: dummy23
  character(len=11) :: spinpar(500)
  character(len=6) :: dummy6
  character(len=4) :: spin(500)
  character(len=3) :: unit
  character(len=1) :: par(500)

  parameter (meV2GHz = 241.798926d0)

  open (7,file=name,status='unknown',form='formatted')
  read(7,'(a24,i3)') dummy24,nvec

  do i=1,3
     read(7,'(a)') string
  end do

  do i=1,nvec
     !read(7,'(a15,d25.10)') level(i),energy(i)
     read(7,'(i4,a11,d25.10)') level(i),spinpar(i),energy(i)
     spin(i) = spinpar(i)(6:9)
     par(i) = spinpar(i)(11:11)
     !write(*,*) spinpar(i)
     !write(*,*) par(i)
     if(spin(i).eq.' 1/2') fspin(i) = 0.5
     if(spin(i).eq.' 3/2') fspin(i) = 1.5
     if(spin(i).eq.' 5/2') fspin(i) = 2.5
     if(spin(i).eq.' 7/2') fspin(i) = 3.5
     if(spin(i).eq.' 9/2') fspin(i) = 4.5
     if(spin(i).eq.'11/2') fspin(i) = 5.5
     if(spin(i).eq.'13/2') fspin(i) = 6.5
     if(spin(i).eq.'15/2') fspin(i) = 7.5
     if(spin(i).eq.'17/2') fspin(i) = 8.5
     if(spin(i).eq.'19/2') fspin(i) = 9.5
     if(spin(i).eq.'21/2') fspin(i) = 10.5

     if(spin(i).eq.'   0') fspin(i) = 0.0
     if(spin(i).eq.'   1') fspin(i) = 1.0
     if(spin(i).eq.'   2') fspin(i) = 2.0
     if(spin(i).eq.'   3') fspin(i) = 3.0
     if(spin(i).eq.'   4') fspin(i) = 4.0
     if(spin(i).eq.'   5') fspin(i) = 5.0
     if(spin(i).eq.'   6') fspin(i) = 6.0
     if(spin(i).eq.'   7') fspin(i) = 7.0
     if(spin(i).eq.'   8') fspin(i) = 8.0
     if(spin(i).eq.'   9') fspin(i) = 9.0
     if(spin(i).eq.'  10') fspin(i) = 10.0
  end do

  do i=1,6
     read(7,'(a)') string
  end do

  do i=1,nvec
     read(7,'(a23,2d20.10,d19.20)') dummy23,nms(i),nmscorr(i),nmsrel(i)
     if(unit.eq.'MHz') then
        nms(i) = nms(i)*1000.d0
        nmscorr(i) = nmscorr(i)*1000.d0
        nmsrel(i) = nmsrel(i)*1000.d0
     else if(unit.eq.'meV') then
        nms(i) = nms(i)/meV2Ghz
        nmscorr(i) = nmscorr(i)/meV2Ghz
        nmsrel(i) = nmsrel(i)/meV2Ghz
     end if
     if(i.lt.nvec) then
        read(7,'(a)') string
        read(7,'(a)') string
        read(7,'(a)') string
     end if
  end do

  do i=1,6
     read(7,'(a)') string
  end do

  do i=1,nvec
     read(7,'(a23,2d20.10,d19.20)') dummy23,sms(i),smscorr(i),smsrel(i)
     !write(*,*) dummy23,sms(i),smscorr(i),smsrel(i)
     if(unit.eq.'GHz') then
        ms(i) = nms(i) + sms(i)
        mscorr(i) = nmscorr(i) + smscorr(i)
        msrel(i) = nmsrel(i) + smsrel(i)
     else if(unit.eq.'MHz') then
        sms(i) = sms(i)*1000.d0
        smscorr(i) = smscorr(i)*1000.d0
        smsrel(i) = smsrel(i)*1000.d0
        ms(i) = nms(i) + sms(i)
        mscorr(i) = nmscorr(i) + smscorr(i)
        msrel(i) = nmsrel(i) + smsrel(i)
     else if(unit.eq.'meV') then
        sms(i) = sms(i)/meV2Ghz
        smscorr(i) = smscorr(i)/meV2Ghz
        smsrel(i) = smsrel(i)/meV2Ghz
        ms(i) = nms(i) + sms(i)
        mscorr(i) = nmscorr(i) + smscorr(i)
        msrel(i) = nmsrel(i) + smsrel(i)
     end if
     if(i.lt.nvec) then
        read(7,'(a)') string
        read(7,'(a)') string
        read(7,'(a)') string
     end if
  end do

  if(eflag.eq.1) then

     do i=1,10+nvec
        read(7,'(a)') string
     end do

     do i=1,nvec
        read(7,'(a23,4d20.10)') dummy23,fs(1,i),fs(2,i),fs(3,i),fs(4,i)
        !write(*,*) dummy23, fs(1,i),fs(2,i),fs(3,i),fs(4,i)
        if(unit.eq.'MHz') then
           fs(1:4,i) = fs(1:4,i)*1000.d0
        else if(unit.eq.'meV') then
           fs(1:4,i) = fs(1:4,i)/meV2Ghz
           !write(*,*) fs(1,i), fs(2,i), fs(3,i), fs(4,i), fs(5,i)
        end if
     end do

     do i=1,5
        read(7,'(a)') string
     end do

     do i=1,nvec
        read(7,'(a23,2d20.10)') dummy23,fs(5,i),fs(6,i)
        if(unit.eq.'MHz') then
           fs(5,i) = fs(5,i)*1000.d0
           fs(6,i) = fs(6,i)*1000.d0
        else if(unit.eq.'meV') then
           fs(5,i) = fs(5,i)/meV2Ghz
           fs(6,i) = fs(6,i)/meV2Ghz
           !write(*,*) fs(1,i), fs(2,i), fs(3,i), fs(4,i), fs(5,i), fs(6,i)
        end if
     end do

  else

     do i=1,6
        read(7,'(a)') string
     end do

     do i=1,nvec
        read(7,'(a23,d20.10)') dummy23,fs(1,i)
        write(*,*) fs(1,i)
        fs(2,i) = 0.0d0
        fs(3,i) = 0.0d0
        fs(4,i) = 0.0d0
        fs(5,i) = 0.0d0
        fs(6,i) = 0.0d0

!        if(unit.eq.'MHz') then
!           fs(1,i) = fs(1,i)*1000.d0
!        else if(unit.eq.'meV') then
!           fs(1,i) = fs(1,i)/meV2Ghz
!           !write(*,*) fs(1,i), fs(2,i), fs(3,i), fs(4,i), fs(5,i)
!        end if
     end do

  end if

  close(7)

  if(rflag.eq.1) then
     nms(:) = nmsrel(:)
     sms(:) = smsrel(:)
     ms(:) = msrel(:)
  end if

  return
end subroutine openfile

SUBROUTINE fermidist(trms,t,w,b20,b40,flag, radmom)
  IMPLICIT NONE
  INTEGER :: m, m2, i, j, k, flag
  double precision :: c0, a, t, w, b20, b40
  double precision :: rmax, h, h2, pi, pi2, rms, trms, drms
  double precision :: dnorm, norm, dinte, inte
  double precision :: dr2, r2, dr4, r4, dr6, r6, dr8, r8, radmom(4)
  double precision :: cv(0:2000), th(0:2000)
  double precision, ALLOCATABLE :: r(:), rho(:)

  parameter (pi = 3.14159265359d0)
  parameter (pi2 = 9.86960440108936d0)

  a = t/(4.d0*log(3.d0))

  rmax = 100.0d0
!  h = 0.1d0
  h = 0.01d0
  m = int(rmax/h)
  m = m + MOD(m,2)

  ALLOCATE ( r(0:m), rho(0:m) )

  c0 = trms
  k = 1
  drms = 1.0d0
  do while (drms >= 1.d-8)
     if(k>1) then
        c0 = trms/rms*c0
     end if
     DO i=0, m
        r(i) = h*i;
        if(abs(b20)<0.001d0.and.abs(b40)<0.001d0) then
           rho(i) = 4.0d0*pi*(1.0d0 + w*r(i)**2.0d0/c0**2.0d0)/(1.0d0+exp((r(i)-c0)/a))
        end if
     END DO

     if(abs(b20) > 0.0d0.or.abs(b40) > 0.0d0) then
        m2 = 360
        h2 = pi/360.d0
        do i=0, m2, 1
           th(i) = i*pi/360.d0
           cv(i) = c0*(1.d0 + b20*sqrt(5.d0/(16.d0*pi))*(3.d0*(cos(th(i)))**2.d0 - 1.d0) &
                + b40*3.d0/(16.d0*sqrt(pi))*(35.d0*(cos(th(i)))**4.d0 - 30.d0*(cos(th(i)))**2.d0 + 3.d0))
        end do

        do i=0, m, 1
           inte = 0.0
           do j=1, m2-1, 2
              dinte = h2/3.0d0 &
                   *((1.0d0 + w*r(i)**2.0d0/c0**2.0d0)/(1.d0 + exp((r(i) - cv(j-1))/a))*sin(th(j-1)) &
                   + 4.0d0*(1.0d0 + w*r(i)**2.0d0/c0**2.0d0)/(1.d0 + exp((r(i) - cv(j))/a))*sin(th(j)) &
                   + (1.0d0 + w*r(i)**2.0d0/c0**2.0d0)/(1.d0 + exp((r(i) - cv(j+1))/a))*sin(th(j+1)))
              inte = inte + dinte
           end do
           rho(i) = 2.0d0*pi*inte
           !write(*,*) r(i), rho(i)
        end do
     end if

     norm = 0.0d0
     r2 = 0.0d0
     r4 = 0.0d0
     r6 = 0.0d0
     r8 = 0.0d0
     DO i=1, m-1, 2
        dnorm = (h/3.0d0)*(rho(i-1)*r(i-1)**2.0d0 + 4.0d0*rho(i)*r(i)**2.0d0 + rho(i+1)*r(i+1)**2.0d0)
        dr2 = (h/3.0d0)*(rho(i-1)*r(i-1)**4.0d0 + 4.0d0*rho(i)*r(i)**4.0d0 + rho(i+1)*r(i+1)**4.0d0)
        dr4 = (h/3.0d0)*(rho(i-1)*r(i-1)**6.0d0 + 4.0d0*rho(i)*r(i)**6.0d0 + rho(i+1)*r(i+1)**6.0d0)
        dr6 = (h/3.0d0)*(rho(i-1)*r(i-1)**8.0d0 + 4.0d0*rho(i)*r(i)**8.0d0 + rho(i+1)*r(i+1)**8.0d0)
        dr8 = (h/3.0d0)*(rho(i-1)*r(i-1)**10.0d0 + 4.0d0*rho(i)*r(i)**10.0d0 + rho(i+1)*r(i+1)**10.0d0)
        norm = norm + dnorm
        r2 = r2 + dr2
        r4 = r4 + dr4
        r6 = r6 + dr6
        r8 = r8 + dr8
     end do
     radmom(1) = r2/norm
     rms = sqrt(radmom(1))
     radmom(2) = r4/norm
     radmom(3) = r6/norm
     radmom(4) = r8/norm
     norm = 1.0d0/norm
     k = k + 1
     if(flag.eq.0) then
        drms = abs(trms-rms)
     else
        drms = 0.0d0
     end if
  end do

  DEALLOCATE ( r, rho )

  return

END SUBROUTINE fermidist

subroutine openlsj(name,nvec,level,spin,par,csflabel,lmax)
  implicit none
  integer :: i, nvec,level(500),lmax
  character(len=4) :: spin(500)
  character(len=1) :: par(500)
  character(len=10) :: spinpar(500)
  character(len=36) :: dummy

  character(len=100) :: name
  character(len=200) :: string
  character(len=200) :: csflabel(500)


  open (7,file=name,status='unknown',form='formatted')

  lmax = 0
  read(7,'(a)') string
  i=1
  do
     read(7,'(i4,a10)',end=99) level(i),spinpar(i)
     if(spinpar(i)(10:10).eq."+".or.spinpar(i)(10:10).eq."-") then
        !read(7,'(i4,a10)') level(i),spinpar(i)
        spin(i) = spinpar(i)(1:4)
        par(i) = spinpar(i)(10:10)
        read(7,'(a36,a)') dummy,csflabel(i)
        if(len(trim(csflabel(i))).gt.lmax) lmax = len(trim(csflabel(i)))
        !write(*,*) level(i),spin(i),' ',par(i),' ',trim(csflabel(i))
        i = i + 1
     !else
     !   spin(i) = spinpar(i)(1:4)
     !   par(i) = spinpar(i)(10:10)
     end if
     !write(*,*) level(i),spin(i),par(i)

     !read(7,'(a36,a)') dummy,csflabel(i)
     !if(len(trim(csflabel(i))).gt.lmax) lmax = len(trim(csflabel(i)))
     !write(*,*) level(i),spin(i),par(i),csflabel(i)

     !read(7,'(a)',end=99) string
     !i = i + 1
  end do
99 continue
  close(7)
  nvec = i-1
  !write(*,*) nvec
end subroutine openlsj
