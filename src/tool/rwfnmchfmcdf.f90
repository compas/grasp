       program ff2gr
!*************************************************************
! Grant and Cowan use 'real' P(r)=r.R(r)
! Froese Fischer and Chernyesheva use:
! P_(r)= P(r)*sqrt(alpha+beta/r)
! FF: alpha=0, beta=1, so P_ff=P(r)/sqrt(r)
! E_gr=0.5*E_cw
! E_ff=    E_cw=EHF()
! E_ch=   -E_cw
!*************************************************************
      PARAMETER (nwf=120,no=220)
      implicit real*8(a-h,o-z)
      DIMENSION pf(no+1),pff(no+1,nwf),rf(no+1),yy2(no+1)
      DIMENSION pg(230),qg(230),rg(230),pgg(230,2)
      DIMENSION e(nwf),az(nwf)
      CHARACTER hfcr*4,line*3,atom*6,term*6
      integer max(nwf),nl(nwf),l(nwf)

!b
!b alpha constant from lib/lib92/setcon.f
!b
!      COMMON/DEF9/CVAC,PI

      write(*,*) 'RWFNMCHFMCDF'
      write(*,*) 'This program converts non-relativistic radial'
      write(*,*) 'orbitals to relativistic ones in GRASP format'
      write(*,*) 'Input file: wfn.inp'
      write(*,*) 'Output file: rwfn.out'


      open (unit=1,file='wfn.inp',status='old',form='UNFORMATTED')
      open (unit=9, file='rwfn.out',status='unknown',form='unformatted')

! u  *******************************************************************
! u  write out initial orbitals Pnl and Qnl as input
! u  for Grant's MCDF program
! u  *******************************************************************

      Z   = 0.d0
      pff = 0.d0 ! JG
      pg  = 0.d0
      qg  = 0.d0

!      do 1 i=1,nwf
!      do 1 j=1,no+1
!         pff(j,i)=0.d0
!1     continue


!      do 2 i=0,230
!         pg(i)=0.d0
!         qg(i)=0.d0
!2     continue


      write(9) 'G92RWF'
!
      do 3 i=1,3
3     line(i:i)=' '

      nomax=0
      do 4 i=1,nwf
         read(1,end=48) atom,term,line,m,z,e(i),ekin,az(i),&
                        (pff(j+1,i),j=1,m)
         max(i)=m+1
         if (m+1.gt. nomax) nomax=m+1
!
! **** read THE NL- and THE L-VALUES
!
         il=1
         call skip (il,line)
         nl(i)= rdnum(il,line)
         l(i)=rdorb(il,line)
4     continue

48    nwf1=i-1

      CALL SETCON

!b    c=137.0359895
!b    halfa=0.5*7.2973531 E-3
      cvac = 137.035999139D0
      halfa=0.5 / cvac

!b end alpha constant

      hfcr='mchf'
      rnt=exp(-4.0625)/z
      h=.0625
      hp=0.
      n=no
      eph = exp (h)
      ett = 1.0

!  *****  construct MCHF grid:
!  *****  GENERATE ARRAYS FOR R WITH A CONSTANT MESH
!  *****  SIZE IN THE LOG(Z*R) VARIABLE

      rho = -4.d0
      hff   = 1.d0/16.d0
      rf(1) = 0.d0
      do 5 i=1,no
         rf(i+1)= exp(rho)/z
         rho = rho + hff
5     continue

!  *****  construct MCDF grid:
!  *****  GENERATE ARRAYS FOR R WITH A CONSTANT MESH
!  *****  SIZE IN THE LOG(Z*R) VARIABLE

      call grid(z,rg)
      h = rg(2)/1000.d0

!*************************************
!** return to original 'real' p(r) ***
!*************************************
      do 6 i=1,nwf
      do 6 j=1,max(i)
         pff(j,i)=dsqrt(rf(j))*pff(j,i)
6     continue

      do 7 m=1,nwf1
         my=max(m)
         do 92 j=1,no+1
            pf(j)=pff(j,m)
92       continue
         if (l(m).eq.0) then
            call spline(rf,pf,221,az(m),0.d0,yy2)
         else
            call spline(rf,pf,221,0.d0,0.d0,yy2)
         endif
         myg=229
         do 94 j=1,230
            if (rg(j).gt.rf(my)) then
               myg = j
               goto 95
            endif
            call splint(rf,pf,yy2,221,rg(j),pg(j))
            call splint(rf,pf,yy2,221,rg(j)-h,pgg(j,1))
            call splint(rf,pf,yy2,221,rg(j)+h,pgg(j,2))
94       continue
95       continue

         if (l(m) .ne. 0) then

!****************************
!** first lower j-value:  ***
!****************************

            pqj=real(l(m)-0.5)
            pqkap=-(pqj+0.5)
            npy=nl(m)
            naky=-nint(pqkap)
            ey=0.5*e(m)

            do 8 j=2,myg
               dpa=(pgg(j,2)-pgg(j,1))/(2.d0*h)
               dpf=(pg(j+1)-pg(j-1))/(rg(j+1)-rg(j-1))
               qg(j)=halfa*(dpa-pqkap*pg(j)/rg(j))
8           continue

            write (9) npy,naky,ey,myg
            write (9) az(m),(pg(i),i=1,myg),(qg(i),i=1,myg)
            write (9) (rg(i),i=1,myg)
         endif

! u  ************************************************
! u       second higher j-value:
! u  ************************************************

         pqj=real(l(m)+0.5)
         pqkap=pqj+0.5
         npy=nl(m)
         naky=-nint(pqkap)
         ey=0.5*e(m)

         do 9 j=2,myg
            dpa=(pgg(j,2)-pgg(j,1))/(2.d0*h)
            dpf=(pg(j+1)-pg(j-1))/(rg(j+1)-rg(j-1))
            qg(j)=halfa*(dpa-pqkap*pg(j)/rg(j))
9        continue

         write (9) npy,naky,ey,myg
         write (9) az(m),(pg(i),i=1,myg),(qg(i),i=1,myg)
         write (9) (rg(i),i=1,myg)

7     continue
      close (1)
      close (9)
      stop

      CONTAINS

      subroutine skip(i,line)
!
! skips irrelevants on a line
!
      character line*3,irrel*8
      data irrel/' ()/\,;:'/

   10 if(index(irrel,line(i:i)) .ne.  0 .and. i .lt. 80) then
        i=i+1
        goto 10
      endif
      end subroutine skip


      function rdnum(I,line)
!
! read AN integer FROM line
!
      integer rdnum
      character*3 line

! ICHAR('I')-ICHAR('0') RETURNS THE integer VALUE OF A character
! ICHAR('0')=48 or 33 DEPendENT ON THE character SET USED, BUT then
! ICHAR('3')=51 or 36, THE FINAL ANSWER BEING CorRECT

   10 if(line(2:2).GE.'0'.and.line(2:2).LE.'?')then
         rdnum=ICHAR(line(2:2))-ICHAR('0')
         I=3
      else
         write(*,*) 'Something wrong with n'
         stop
      endif
      end function

      function rdorb(i,line)

!
! converts character to orbital angular momentum
!
      integer rdorb
      character lorbu*1,lorbl*1,line*3

      dimension lorbu(0:10),lorbl(0:10)
      data lorbu/'S','P','D','F','G','H','I','K','L','M','N'/
      data lorbl/'s','p','d','f','g','h','i','k','l','m','n'/

      rdorb=0
      do 1 ii=0,10
1     if (line (i:i) .eq. lorbl(ii) .or.      &
          line (i:i) .eq. lorbu(ii)) rdorb=ii
      I=I+1
      end  function


      subroutine GRID(z,rr)
!
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RR(230)

      RNT = EXP (-65.0D00/16.0D00) / Z
      H = 0.5D00**4
      N = MIN (220,400)
      NP10 = N+10
      RR(1) = 0.0D00
      EPH = EXP (H)
      ETT = 1.0D00
!
!   Set up the arrays R, RP, RPOR
!
      DO 1 I = 2,NP10
         ETT = EPH*ETT
         ETTM1 = ETT-1.0D00
         RR(I) = RNT*ETTM1
 1    CONTINUE
      return
      end SUBROUTINE

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+ 1)                        &
       -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*   &
        u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END SUBROUTINE

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) print *, ' bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
        (h**2)/6.d0
      return
      END  SUBROUTINE
      end Program
