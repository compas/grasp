      PROGRAM WFNPLOT
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER AT*6,TT*6,EL1*3,EL*3,NEW*3,INPUT*24,OUTPUT*24,  &
          ATOM*6,TERM*6,XA*3,OUTPUT2*24,leg*50,                &
          format_string*50,name*100
      DIMENSION P(220,30),EL(30),R(220),R2(220)
!
      write(*,*) '****************************************************'
      write(*,*) 'Program wfnplot writes HF/MCHF radial wave functions'
      write(*,*) 'to following output files:'
      write(*,*)
      write(*,*) 'Matlab/GNU Octave file "octave_name.m"'
      write(*,*) 'Xmgrace file "xmgrace_name.agr"'
      write(*,*)
      write(*,*) 'Input file:  name.w'
      write(*,*)
      write(*,*) 'To plot orbital: press enter'
      write(*,*) 'To remove orbital: type "d" or "D" and press enter'
      write(*,*)
      write(*,*) '                             Jorgen Ekman Jun 2015'
      write(*,*) '****************************************************'

      inc = 1

      WRITE(*,*) 'Name of state:'
      READ(*,*) name
      OPEN(3,FILE=trim(name)//'.w',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(4,FILE='octave_'//trim(name)//'.m',STATUS='UNKNOWN')
      OPEN(8,FILE='xmgrace_'//trim(name)//'.agr',STATUS='UNKNOWN')
      write(*,*)
      write(*,*) 'To have r on x-axis: type "y" otherwise "n" for sqrt(r '
      write(*,*)
      read(*,*) XA

      IUF=3
      nwf=1
      MM = 0
2     P(1,nwf) = 0.D0
      READ(IUF,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(P(J,NWF),J=2,M+1)
      WRITE(6,'(2x,A,A,$)') EL1,' = '
      READ(5,'(A)') NEW
      IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
         IF ( NEW .NE. '   ') THEN
            EL1 = NEW
         ENDIF
         EL(NWF) = EL1
	 NWF = NWF + 1
         Z = ZT
	 ATOM = AT
	 TERM = TT
	 MM = MAX(M,MM)
      END IF
      GO TO 2
5     CLOSE(UNIT=3)
      RHO = -4.0
      H = 1/16.
      R(1) = 0.d0
      R2(1) = 0.d0
      DO 10 J = 2,220
         R(J) = EXP(RHO)/Z
         R2(J) = SQRT(R(J))
	 RHO = RHO + H
10    CONTINUE
!
!	LIST TABLES OF RADIAL FUNCTIONS
!
      NWF = NWF -1
      write(4,*) '% sqrt(r)    P(nl;r)'
      write(4,*) 'clf'

      if(XA.eq.'y  ') then
         write(8,*) '@    xaxis  label "r"'
      else
         write(8,*) '@    xaxis  label "sqrt(r)"'
      end if
      write(8,*) '@    yaxis  label "P(r)"'
      write(8,*) '# sqrt(r)    P(nl;r)'

      DO 20 I = 1,NWF

         if(I-1.lt.10) then
            format_string = "(A,I1,3A)"
         else
            format_string = "(A,I2,3A)"
         end if
         write (leg,format_string)  &
             '@    s',I-1,'  legend  "',EL(I),'"'
         write(8,*) trim(leg)
         write(8,*) '# ',EL(I)

	 write(4,*) 'P = ['
	 DO 21 J = 1,MM,inc
	    Pwave = P(J,I)*R2(J)
            IF(XA.EQ.'y  ') THEN
               IF (abs(Pwave) .gt. 0.0005 .OR. J.EQ.1) 	&
                   WRITE(4,'(F10.4,F12.3)') R(J),Pwave
               WRITE(8,'(F10.4,F12.3)') R(J),Pwave
            ELSE
               IF (abs(Pwave) .gt. 0.0005 .OR. J.EQ.1) 	&
                   WRITE(4,'(F10.4,F12.3)') R2(J),Pwave
               WRITE(8,'(F10.4,F12.3)') R2(J),Pwave
            END IF
 21      CONTINUE
         write(8,*)
	 write(4,*) '];'
         write(4,*) 'plot(P(:,1), P(:,2))'
         write(4,*) 'hold all'
 20   CONTINUE
      IF(XA.EQ.'y  ') THEN
         write(4,*) 'xlabel (''r'', ''fontsize'', 12)'
      ELSE
         write(4,*) 'xlabel (''sqrt(r)'', ''fontsize'', 12)'
      END IF
      write(4,*) 'ylabel (''P(r)'', ''fontsize'', 12)'
      write(4,*) 'grid on'
      write(4,101,advance='no') 'legend('
      DO 22 I = 1,NWF
         IF(I.LT.NWF) THEN
            write(4,102,advance='no') '''',EL(I),''','
         ELSE
            write(4,102,advance='no') '''',EL(I),''')'
         END IF
 22   CONTINUE
 101  format(a)
 102  format(3a)
      END
