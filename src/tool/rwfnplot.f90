      PROGRAM rwfnplot
      IMPLICIT NONE

      INTEGER, PARAMETER:: NPTS0=5000

      DOUBLE PRECISION pg(NPTS0), qg(NPTS0), rg(NPTS0), pgg(NPTS0)
      DOUBLE PRECISION rg2(NPTS0)
      DOUBLE PRECISION energy, a0
      CHARACTER        title*6, orbl*4, nnstr*2, new*3, xa*3, el*5
      CHARACTER        leg*50, format_string*50, name*100
      DIMENSION        el(30)
      INTEGER          np, lp, jp, nn, laky, ll, jj, npts, j, i, nwf

      write(*,*) 'RWFNPLOT'
      write(*,*) 'Program to generate Matlab/GNU Octave and'
      write(*,*) 'Xmgrace files that plot radial orbitals'
      write(*,*) 'Input file:  name.w'
      write(*,*) 'Output files: octave_name.m, xmgrace_name.agr'
      write(*,*)
      write(*,*) 'To plot orbital: press enter'
      write(*,*) 'To remove orbital: type "d" or "D" and press enter'
      write(*,*)
      write(*,*) '                             Jorgen Ekman Jun 2015'
      write(*,*)

      WRITE(*,*) 'Name of state:'
      READ(*,*) name
      OPEN(3,FILE=trim(name)//'.w',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(4,FILE='octave_'//trim(name)//'.m',STATUS='UNKNOWN')
      OPEN(8,FILE='xmgrace_'//trim(name)//'.agr',STATUS='UNKNOWN')

      write(*,*)
      write(*,*) 'To have r on x-axis: type "y" otherwise "n" for sqrt(r)'
      read(*,*) xa
      write(*,*)
      i = 1

      write(4,*) '%     r                   P(r)                   Q(r)'

! Xmgrace stuff
      if(xa.eq.'y  ') then
         write(8,*) '@    xaxis  label "r"'
      else
         write(8,*) '@    xaxis  label "sqrt(r)"'
      end if
      write(8,*) '@    yaxis  label "P(r)"'
      write(8,*) '#     r                   P(r)                   Q(r)'

      write(4,*) 'clf'
      READ(3) title
      IF (title .NE. 'G92RWF') THEN   ! Extra safety
         PRINT *, 'title = ', title, 'does not match G92RWF'
         STOP
      ENDIF

      DO
         READ(3, END = 20) nn, laky, energy, npts
         write(nnstr,'(I2)') nn
         IF (laky .GT. 0) THEN
            ll = laky
            jj = -1
         ELSEIF (laky .LE. -1) THEN
            ll = -laky - 1
            jj = 1
         ELSE
            WRITE(*,*)'Unexpected case in reading mcdf.w'
            STOP
         ENDIF

         if(jj.eq.1) then
            if(ll.eq.0)  orbl = nnstr//'s '
            if(ll.eq.1)  orbl = nnstr//'p '
            if(ll.eq.2)  orbl = nnstr//'d '
            if(ll.eq.3)  orbl = nnstr//'f '
            if(ll.eq.4)  orbl = nnstr//'g '
            if(ll.eq.5)  orbl = nnstr//'h '
            if(ll.eq.6)  orbl = nnstr//'i '
            if(ll.eq.7)  orbl = nnstr//'j '
            if(ll.eq.8)  orbl = nnstr//'k '
            if(ll.eq.9)  orbl = nnstr//'l '
         else
            if(ll.eq.1)  orbl = nnstr//'p-'
            if(ll.eq.2)  orbl = nnstr//'d-'
            if(ll.eq.3)  orbl = nnstr//'f-'
            if(ll.eq.4)  orbl = nnstr//'g-'
            if(ll.eq.5)  orbl = nnstr//'h-'
            if(ll.eq.6)  orbl = nnstr//'i-'
            if(ll.eq.7)  orbl = nnstr//'j-'
            if(ll.eq.8)  orbl = nnstr//'k-'
            if(ll.eq.9)  orbl = nnstr//'l-'
         end if

         IF (npts .GT. NPTS0) THEN
            WRITE(*,*) 'df2hf: npts .GT. NPTS0'
            STOP
         ENDIF

         READ(3) a0, (pg(j), j=1,npts), (qg(j), j=1,npts)
         READ(3) (rg(j), j=1,npts)

         WRITE(*,'(2x,A,A,$)') orbl,' = '
         READ(*,'(A)') new
         IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
            el(i) = orbl

            ! Xmgrace stuff
            if(i-1.lt.10) then
               format_string = "(A,I1,3A)"
            else
               format_string = "(A,I2,3A)"
            end if
            write (leg,format_string)   &
                 '@    s',i-1,'  legend  "',el(i),'"'
            write(8,*) trim(leg)
            write(8,*) '# ',el(i)

            i = i + 1
            write(4,*) 'P = ['
            DO j = 1, npts
               rg2(j) = sqrt(rg(j))
               if(xa.eq.'y  ') then
                  if (abs(pg(j)) .gt. 0.0005 .OR. j.EQ.1) 	&
                       WRITE (4, '(3D20.10)') rg(j), pg(j), qg(j)
                  WRITE (8, '(3D20.10)') rg(j), pg(j), qg(j)
               else
                  if (abs(pg(j)) .gt. 0.0005 .OR. j.EQ.1) 	 &
                  WRITE (4, '(3D20.10)') rg2(j), pg(j), qg(j)
                  WRITE (8, '(3D20.10)') rg2(j), pg(j), qg(j)
               end if
            ENDDO
            write(8,*)
            write(4,*) '];'
            write(4,*) 'plot(P(:,1), P(:,2))'
            write(4,*) 'hold all'
         ENDIF
      ENDDO
 20   CONTINUE
      nwf = i - 1
      if(xa.eq.'y  ') then
!         write(4,*) 'xlabel (''r'', ''fontsize'', 12)'
         write(4,*) 'xlabel (''r'')'
      else
!         write(4,*) 'xlabel (''sqrt(r)'', ''fontsize'', 12)'
         write(4,*) 'xlabel (''sqrt(r)'')'
      end if
!      write(4,*) 'ylabel (''P(r)'', ''fontsize'', 12)'
      write(4,*) 'ylabel (''P(r)'')'
      write(4,*) 'grid on'
      write(4,101,advance='no') 'legend('
      DO i = 1,nwf
         IF(i.LT.nwf) THEN
            write(4,102,advance='no') '''',el(i),''','
         ELSE
            write(4,102,advance='no') '''',el(i),''')'
         END IF
      ENDDO
 101  format(a)
 102  format(3a)
      PRINT *, ' FINISHED .....'
      STOP
      END
