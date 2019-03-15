!***********************************************************************
!                                                                      *
      PROGRAM TRANSTABLE
!                                                                      *
!   This program reads the output from biotra2 together with spectral  *
!   designation and energies to produce a latex table                  *
!                                                                      *
!   Written by Per Jonsson, Malmo University, November 2010            *
!                                                                      *
!***********************************************************************
!
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      PARAMETER (NTRANS = 100000)
      INTEGER ILAB(NTRANS),GFPRINT,IU(NTRANS),IL(NTRANS),INDX(NTRANS)
      INTEGER ASCII
      CHARACTER*80 ROW,FILENAME1,FILENAME2,LATEX(NTRANS),ULATEX(NTRANS)
      CHARACTER*80 LLATEX(NTRANS)
      CHARACTER*80 LATEXFILE1,LATEXFILE2
      CHARACTER*2 FU(NTRANS),FL(NTRANS),GAUGE(NTRANS),FLAB(NTRANS)
      CHARACTER*2 EM(NTRANS)
      CHARACTER*4 JU(NTRANS),JL(NTRANS),JLAB(NTRANS),PU(NTRANS)
      CHARACTER*4 PL(NTRANS)
      CHARACTER*1 ANS,PAR
      CHARACTER*14 FILE(NTRANS),FILE1,FILE2,FILEU(NTRANS),FILEL(NTRANS)
      CHARACTER*14 DUMMY
      CHARACTER*26 ROW1
      CHARACTER*100 STRING,LC
      CHARACTER*3 P(NTRANS)
      REAL*8 UENERGY(NTRANS),LENERGY(NTRANS),ENERGY(NTRANS)
      REAL*8 S(NTRANS),GF(NTRANS),DELTAE(NTRANS),R(NTRANS)
      REAL*8 DELTAA(NTRANS),AL(NTRANS),AV(NTRANS),TL(NTRANS),TV(NTRANS)
      REAL*8 DTMEAN,FRAC
!

      WRITE(*,*)
      WRITE(*,*) ' RTABTRANS2    '
      WRITE(*,*) ' This program reads energy label data and transition'
      WRITE(*,*) ' data and creates transition and lifetime tables in '
      WRITE(*,*) ' LaTeX or ASCII format. An Octave file with a       '
      WRITE(*,*) ' scatterplot of dT and 10log(A) is also produced    '
      WRITE(*,*)
      WRITE(*,*) ' Energy label data are given in the file energylabel'
      WRITE(*,*) ' created by the rtabtrans1 program                  '
      WRITE(*,*) ' Transition data file can be conctenated *.t or *.ct'
      WRITE(*,*) ' files. '
      WRITE(*,*)
      WRITE(*,*) ' Input files: energylabel.latex(ascii),             '
      WRITE(*,*) ' transitiondatafile                                 '
      WRITE(*,*) ' Output files: transitiontable.tex(txt),            '
      WRITE(*,*) '               lifetimetable.tex(txt),              '
      WRITE(*,*) '               scatterplot.m                        '
      WRITE(*,*)
      WRITE(*,*) ' Give the name of the transition data file          '
      READ(*,*) FILENAME1

      WRITE(*,*) ' Energy label file in LaTeX or ASCII format (0/1)?  '
      READ(*,*) ASCII
      IF (ASCII.EQ.1) THEN
         WRITE(*,*) ' Left column string to denote the ion, e.g. Fe X'
         READ(*,'(A)') LC
      END IF

      IF (ASCII.EQ.1) THEN
         FILENAME2 = 'energylabel.ascii'
         LATEXFILE1 = 'transitiontable.txt'
         LATEXFILE2 = 'lifetimetable.txt'
      ELSE
         FILENAME2 = 'energylabel.latex'
         LATEXFILE1 = 'transitiontable.tex'
         LATEXFILE2 = 'lifetimetable.tex'
      END IF

      WRITE(*,*) ' Give cut-off for printing A values '
   	READ(*,*) CUTOFF
      WRITE(*,*) ' Give fraction of accumulated A value for upper level'
      write(*,*) ' for printing A value of a transition '
      READ(*,*) FRAC
      WRITE(*,*) ' Transition data wavelength sorted? '
      READ(*,'(A)') ANS
      IF ((ANS.EQ.'y').OR.(ANS.EQ.'Y')) THEN
         NSORT = 1
      ELSE
         NSORT = 0
      END IF
      WRITE(*,*) ' Give number of decimals for wavelength (1,...6) '
      READ(*,*)  IDIGITS
      WRITE(*,*)

      OPEN(36,FILE=FILENAME1,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(37,FILE=FILENAME2,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(38,FILE=LATEXFILE1,FORM='FORMATTED',STATUS = 'UNKNOWN')
      OPEN(39,FILE=LATEXFILE2,FORM='FORMATTED',STATUS = 'UNKNOWN')
      OPEN(40,FILE='scatterplot.m',FORM='FORMATTED',STATUS = 'UNKNOWN')

!---- Read energylabelfile --------------------

       NFOUND = 0
       DO
          READ(37,'(A)') STRING
          IF (STRING(1:5).EQ.'-----') THEN
             NFOUND = NFOUND + 1
          END IF
          IF (NFOUND.EQ.3) EXIT
       END DO


      I = 1
      NLEN = 0
      DO
        READ(37,200,END=91) NO,ILAB(I),JLAB(I),PAR,ENERGY(I), &
              DENER,FILE(I),LATEX(I)
        P(I) = '   '
        P(I)(2:2) = PAR
        IF (LEN_TRIM(LATEX(I)).GT.NLEN) NLEN = LEN_TRIM(LATEX(I))
        JLAB(I) = ADJUSTR(JLAB(I))
        I = I + 1
      END DO

   91 CONTINUE
      N = I - 1

!---- Check that there are no two quantum labels that are the same

      NQEQUAL = 0
      DO I = 1,N
         DO K = I+1,N
            IF ( (LATEX(I).EQ.LATEX(K)).AND.(JLAB(I).EQ.JLAB(K)) ) THEN
               WRITE(*,*) 'Quantum labels for states',I,K,'equal'
               NQEQUAL = NQEQUAL + 1
            END IF
         END DO
      END DO
      IF (NQEQUAl.GT.0) THEN
         WRITE(*,*) 'Some quantum states have equal labels'
         WRITE(*,*) 'Update energylabel file and rerun'
         STOP
      END IF


!---- Read transition files ---------------------------------------

      K = 1
      DO
   98   CONTINUE

!---- Start reading the file. Find out if transition between configurations
!     or within configuration

        READ(36,600) ROW1
        IF (LEN(TRIM(ROW1)).GT.22) THEN
          N12 = 2
          READ(36,500) FILE1
          READ(36,500) FILE2
    	  ELSE
          N12 = 1
          READ(36,501) FILE1
        END IF

!---- Loop over multipolarities in the file

        DO
  101     CONTINUE

!---- Find out if electric or magnetic case

          DO I = 1,3
            READ(36,'(A)') ROW
          END DO

          IF (ROW(2:9).EQ.'Electric') THEN
            NEM = 1
          ELSE
            NEM = 2
          END IF
   	    IF (ROW(16:16).EQ.'1') THEN
	         NEMO = 1
	       ELSE
	         NEMO = 2
	       END IF

!---- Continue to read until the rate information comes

          DO I = 1,4
            READ(36,'(A)') ROW
          END DO

!---- Start reading the rates for the found polarity (electric or magnetic)
          DO
            READ(36,300) FU(K),IU(K),JU(K),PU(K),FL(K),IL(K),JL(K), &
      	                PL(K),DELTAE(K),GAUGE(K),AV(K),GF(K),S(K)

!---- If magnetic then save rate both in AV(K) and AL(K)
!     If electric AL(K) will be overwritten as the next line is read

            AL(K) = AV(K)
            R(K) = 0.D0

            IF (NEM.EQ.1) THEN
              READ(36,301) GAUGE(K),AL(K),GF(K),S(K)
              R(K) = dabs(AL(K)-AV(K))/maxval([AL(K),AV(K)])
            END IF
            IF (FU(K).EQ.'f1') THEN
              FILEU(K) = FILE1
              FILEL(K) = FILE2
	    ELSE IF (FU(K).EQ.'f2') THEN
	      FILEU(K) = FILE2
              FILEL(K) = FILE1
	    ELSE
	      FILEU(K) = FILE1
              FILEL(K) = FILE1
	    END IF

!-----Set polarity for case

            IF ((NEM.EQ.1).AND.(NEMO.EQ.1)) THEN
		        EM(K) = 'E1'
	         ELSE IF ((NEM.EQ.1).AND.(NEMO.EQ.2)) THEN
              EM(K) = 'E2'
            ELSE IF ((NEM.EQ.2).AND.(NEMO.EQ.1)) THEN
              EM(K) = 'M1'
	         ELSE
	           EM(K) = 'M2'
	         END IF

!---- Try and read another line
!     Four cases occur: (1) end of file in which case we close the preset file
!                       (2) we can read a line and continue reading in the normal way
!                       (3) we read a blank line which indicates that we have one
!                            more polarity to read
!                       (4) we can read a file starting with 'Transition' in which case
!                           we have more to read

            READ(36,700,END=99) ROW
            K = K + 1
            IF (K.EQ.NTRANS) THEN
               WRITE(*,*) 'Too many transitions'
               WRITE(*,*) 'Increase NTRANS and recompile'
               STOP
            END IF
            BACKSPACE 36

!---- More polarity

            IF (LEN(TRIM(ROW)).EQ.0) THEN
			     GOTO 101
	         END IF

!---- More transitions but within new files

            IF (LEN(TRIM(ROW)).LT.30) THEN
			     GOTO 98
	         END IF

!---- If not the above cases just continue with another row in the normal way

          END DO
	     END DO
      END DO

   99 CONTINUE


!---- Write header to LaTeX or ASCII  lifetime table ----------------

      IF (ASCII.EQ.1) THEN
         WRITE(39,'(A)') 'Ion  J p  state     Tl(s)      Tv(s)'
      ELSE
         write(39,'(A)') '\documentclass[10pt]{article}'
         write(39,'(A)') '\usepackage{longtable}'
         write(39,'(A)') '\begin{document}'
	      WRITE(39,'(A)') '\begin{longtable}{lll} \hline'
         WRITE(39,'(A)') 'State & $\tau_l$ & $\tau_v$  \\ \hline'
      END IF

!---- Compute and print lifetimes for all the states ---------------

      DO J = 1,N
        TL(J) = 0.D0
        TV(J) = 0.D0
        DO I = 1,N
          DO L = 1,K
            IF ((FL(L) == FL(L)).AND.(IL(L) == ILAB(I)).AND.          &
                (JL(L) == JLAB(I)).AND.(FILEL(L) == FILE(I)).AND.     &
                (FU(L) == FU(L)).AND.(IU(L) == ILAB(J)).AND.          &
                (JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
              TL(J) = TL(J) + AL(L)
              TV(J) = TV(J) + AV(L)
            END IF
          END DO
        END DO
        IF (TL(J).GT.0.D0) THEN
          IF (ASCII.EQ.1) THEN
            WRITE(39,520) TRIM(LC),'  ',JLAB(J),P(J),LATEX(J)(1:NLEN),&
               '   ',1.D0/TL(J),'   ',1.D0/TV(J)
          ELSE
            WRITE(39,420) LATEX(J)(1:NLEN),' & ',1.D0/TL(J),          &
                                         ' & ',1.D0/TV(J),'\\'
          END IF
        END IF
      END DO

      IF (ASCII.EQ.0) THEN
         write(39,'(A)') '\hline\\'
         write(39,'(A)') '\caption{Lifetimes in s.}'
         write(39,'(A)') '\end{longtable}'
         write(39,'(A)') '\end{document}'
      END IF

!---- Write header to LaTeX or ASCII table ----------------------

      IF (ASCII.EQ.1) THEN
         WRITE(38,'(A)') 'Ion Ju pu ustate Jl pl lstate EM Delta E(cm-1) &
         lambda(AA) A(s-1) gf  dT'
      ELSE
         write(38,'(A)') '\documentclass[10pt]{article}'
         write(38,'(A)') '\usepackage{longtable}'
         write(38,'(A)') '\begin{document}'
	      WRITE(38,'(A)') '\begin{longtable}{lllrllll} \hline'
	      WRITE(38,'(A)') 'Upper & Lower & EM & ',                   &
         '$\Delta E$ (cm$^{-1}$) & $\lambda$ (\AA) & ',                  &
         '$A$ (s$^{-1}$) & $gf$ & $dT$                                   &
               \\ \hline'
      END IF

!---- Write header to scatterplot file -------------------------

      write(40,'(A)') 'A = ['

!---- Sort energies --------------------------------------------

      INDX = 0

      DO J = 1,N
        DO I = 1,N
          DO L = 1,K
            IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.         &
               (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J))            &
                .AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
                DELTAA(L) = 1.0E+8/DELTAE(L)
	         END IF
          END DO
        END DO
      END DO

      call HPSORT(K,DELTAA,indx)

      DTMEAN = 0.D0
      DTCOUNT = 0.D0

      IF (NSORT.EQ.0) THEN

!---- Write transition data to LaTeX file. ----------------------
!     Order according to the order in energy file

        DO J = 1,N
          DO I = 1,N
            DO L = 1,K
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.        &
                (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J))            &
                .AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
	             IF ((AL(L).GT.CUTOFF).OR.(AL(L).GT.FRAC*TL(J))) THEN
!               IF (AL(L).GT.CUTOFF) THEN
                DELTAA(L) = 1.0E+8/DELTAE(L)
                IF (EM(L) .EQ. 'M1' .OR.  EM(L) .EQ. 'M2') THEN
                 IF (ASCII.EQ.1) THEN
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,526)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,525)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,524)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,523)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,522)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,521)TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 ELSE
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,426)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,425)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,424)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,423)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,422)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,421)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN) &
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',              &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 END IF
                ELSE
                 WRITE(40,*) AL(L), R(L)
                 DTMEAN = DTMEAN + R(L)
                 DTCOUNT = DTCOUNT + 1.D0
                 IF (ASCII.EQ.1) THEN
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,516) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,515) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,514) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,513) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,512) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,511) TRIM(LC),'  ',JLAB(J),P(J),   &
                     LATEX(J)(1:NLEN),'    ',                   &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),             &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',      &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 ELSE
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,416)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,415)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,414)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,413)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,412)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,411)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 END IF
                END IF
	       END IF
	      END IF
            END DO
          END DO
        END DO

      ELSE

!---- Write transition data to LaTeX file. ----------------------
!     Order according to wave length

        DO LL = 1,K
          L = INDX(LL)
          DO J = 1,N
            DO I = 1,N
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.    &
                (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J))        &
                .AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
	       IF ((AL(L).GT.CUTOFF).OR.(AL(L).GT.FRAC*TL(J))) THEN
                  DELTAA(L) = 1.0E+8/DELTAE(L)
                IF (EM(L) .EQ. 'M1' .OR.  EM(L) .EQ. 'M2') THEN
                 IF (ASCII.EQ.1) THEN
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,526)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,525)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,524)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,523)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,522)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,521)TRIM(LC),'  ',JLAB(J),P(J), &
                     LATEX(J)(1:NLEN),'    ',                &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),          &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',   &
                    DELTAA(L),'  ',AL(L),'  ',GF(L)
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 ELSE
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,426)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,425)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,424)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,423)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,422)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,421)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &','\\'
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 END IF
                ELSE
                 WRITE(40,*) AL(L), R(L)
                 DTMEAN = DTMEAN + R(L)
                 DTCOUNT = DTCOUNT + 1.D0
                 IF (ASCII.EQ.1) THEN
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,516) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,515) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,514) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,513) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,512) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,511) TRIM(LC),'  ',JLAB(J),P(J),  &
                     LATEX(J)(1:NLEN),'    ',                  &
                     JLAB(I),P(I),LATEX(I)(1:NLEN),            &
                    '   ',EM(L),'   ',INT(DELTAE(L)),'  ',     &
                    DELTAA(L),'  ',AL(L),'  ',GF(L),'  ',0.1*R(L)
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 ELSE
                  IF(IDIGITS .EQ. 6) THEN
                    WRITE(38,416)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 5) THEN
                    WRITE(38,415)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 4) THEN
                    WRITE(38,414)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 3) THEN
                    WRITE(38,413)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 2) THEN
                    WRITE(38,412)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSEIF (IDIGITS .EQ. 1) THEN
                    WRITE(38,411)LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN)&
                    ,' & ',EM(L),' & ',INT(DELTAE(L)),' &',             &
                    DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
                  ELSE
                    WRITE(*,*) ' Number of decimals must be 1,...6 '
                    STOP
                  END IF
                 END IF
	             END IF
	            END IF
	           END IF
            END DO
          END DO
        END DO

      END IF

      WRITE(*,*)
      WRITE(*,*) ' Mean dT',DTMEAN/DTCOUNT

      IF (ASCII.EQ.0) THEN
         write(38,'(A)') '\hline\\'
         write(38,'(A)') '\caption{Transition data}'
         write(38,'(A)') '\end{longtable}'
         write(38,'(A)') '\end{document}'
      END IF

      write(40,'(A)') '];'
      write(40,'(A)') "semilogx(A(:,1),A(:,2),'+')"
      write(40,'(A)') "title('Scatterplot of dT and A (s^{-1})')"
      write(40,'(A)') "xlabel('A (s^{-1})')"
      write(40,'(A)') "ylabel('dT')"

      WRITE(*,*)
      WRITE(*,*) ' Program finished. The transition tables in latex'
      WRITE(*,*) ' have been written to file '

 200  FORMAT(2I3,1X,A4,1x,A1,2X,F14.7,F12.2,2X,A14,6X,A)
!cjb  format for highly charged ions in rtransition90/printa.f90 changed
!cjb  format for highly charged ions in rtransition90/printa.f90 F11.2 -> F13.2
!cjb  format for highly charged ions here 1D12.6 -> 1D14.6
!300  FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,1D12.6,A2,1P,   &
 300  FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,1D14.6,A2,1P,   &
           D13.5,2D12.4)
 301  FORMAT(43X,A2,1P,D13.5,2D12.4)
 416  FORMAT(6A,I11,A,F20.6,A,1P,E11.3,A,E11.3,A,F11.3,A)
 415  FORMAT(6A,I11,A,F20.5,A,1P,E11.3,A,E11.3,A,F11.3,A)
 414  FORMAT(6A,I11,A,F20.4,A,1P,E11.3,A,E11.3,A,F11.3,A)
 413  FORMAT(6A,I11,A,F20.3,A,1P,E11.3,A,E11.3,A,F11.3,A)
 412  FORMAT(6A,I11,A,F20.2,A,1P,E11.3,A,E11.3,A,F11.3,A)
 411  FORMAT(6A,I11,A,F20.1,A,1P,E11.3,A,E11.3,A,F11.3,A)
 426  FORMAT(6A,I11,A,F20.6,A,1P,E11.3,A,E11.3,A,A)
 425  FORMAT(6A,I11,A,F20.5,A,1P,E11.3,A,E11.3,A,A)
 424  FORMAT(6A,I11,A,F20.4,A,1P,E11.3,A,E11.3,A,A)
 423  FORMAT(6A,I11,A,F20.3,A,1P,E11.3,A,E11.3,A,A)
 422  FORMAT(6A,I11,A,F20.2,A,1P,E11.3,A,E11.3,A,A)
 421  FORMAT(6A,I11,A,F20.1,A,1P,E11.3,A,E11.3,A,A)
 420  FORMAT(2A,1P,E11.4,A,E11.4,A)
 520  FORMAT(6A,1P,E11.4,A,E11.4)
 500  FORMAT(6X,A14)
 501  FORMAT(5X,A14)
 600  FORMAT(A26)
 700  FORMAT(A80)
 516  FORMAT(12A,I11,A,F20.6,A,1P,E11.3,A,E11.3,A,F11.3)
 515  FORMAT(12A,I11,A,F20.5,A,1P,E11.3,A,E11.3,A,F11.3)
 514  FORMAT(12A,I11,A,F20.4,A,1P,E11.3,A,E11.3,A,F11.3)
 513  FORMAT(12A,I11,A,F20.3,A,1P,E11.3,A,E11.3,A,F11.3)
 512  FORMAT(12A,I11,A,F20.2,A,1P,E11.3,A,E11.3,A,F11.3)
 511  FORMAT(12A,I11,A,F20.1,A,1P,E11.3,A,E11.3,A,F11.3)
 526  FORMAT(12A,I11,A,F20.6,A,1P,E11.3,A,E11.3)
 525  FORMAT(12A,I11,A,F20.5,A,1P,E11.3,A,E11.3)
 524  FORMAT(12A,I11,A,F20.4,A,1P,E11.3,A,E11.3)
 523  FORMAT(12A,I11,A,F20.3,A,1P,E11.3,A,E11.3)
 522  FORMAT(12A,I11,A,F20.2,A,1P,E11.3,A,E11.3)
 521  FORMAT(12A,I11,A,F20.1,A,1P,E11.3,A,E11.3)

      CONTAINS

      SUBROUTINE HPSORT(N,RA,IND)
      double precision RA(N)
      integer IND(N)
      L=N/2+1
      IR=N
      do I=1,N
         IND(I) = I
      end do

  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10    continue
      if(L > 1)then
         L=L-1
         RRA=RA(L)
         IRRA = IND(L)  ! je
      else
         RRA=RA(IR)
         IRRA = IND(IR) ! je
         RA(IR)=RA(1)
         IND(IR)=IND(1)  !je
         IR=IR-1
         if(IR.eq.1)then
            RA(1)=RRA
            IND(1)=IRRA !je
            return
         end if
      end if
      I=L
      J=L+L
20    if(J.le.IR)then
         if(J < IR)then
            if(RA(J) < RA(J+1))  J=J+1
         end if
         if(RRA < RA(J))then
            RA(I)=RA(J)
            IND(I)=IND(J) !je
            I=J; J=J+J
         else
            J=IR+1
         end if

         goto 20
      end if
      RA(I)=RRA
      IND(I)=IRRA
      goto 10
      END SUBROUTINE
      END PROGRAM TRANSTABLE
