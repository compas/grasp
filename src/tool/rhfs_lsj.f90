!************************************************************************
!                                                                      *
      PROGRAM rhfs_lsj
!                                                                      *
!     Program for displaying atomic hyperfine structure parameters     *
!     and Lande g-factors. The program reads the LSJ classification    *
!     file for labeling purposes.
!                                                                      *
!     Per Jonsson and Gediminas Gaigalas                               *    
!                                               August 2011            *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER:: JMax = 22      ! max J value, see JFraction !
      INTEGER, PARAMETER:: ndim = 20000   ! max number of states
      INTEGER, PARAMETER:: maxFile = 1000 ! max number of files 
!
      CHARACTER(LEN=80) strInFile(maxFile), strFile
      CHARACTER*1 iaspa(ndim), PlusMinus(-1:1),ans    ! Parity
      CHARACTER*1 Lev_par(ndim),SWAPP                 ! Parity
      CHARACTER*4 iatjp(ndim), JFraction(1:2*Jmax+1)  ! J
      CHARACTER*4 Lev_J(ndim),SWAPJ                   ! J
      CHARACTER g92mix*6
      CHARACTER*64 string_CSF(ndim), string_PRN(ndim),SWAPCSF ! String in LSJ
      CHARACTER*256 util_lbl_file, hfs_file, output_file
      CHARACTER*32 DUMMY0
      CHARACTER*64 DUMMY
!
      INTEGER i, j, iargc, ios, ncountState, nFile, mFile
      INTEGER nelec, ncftot, nw, nvectot, nvecsiz, nblock, jblock
      INTEGER nb, ncfblk, nevblk, iiatjp, iiaspa
      INTEGER K, Iprint, IMaxCount, nsort
      INTEGER ivec(ndim), indx(ndim), Lev_POS(ndim), izero
      INTEGER MAX_STRING_LENGTH,N
!
      DOUBLE PRECISION eav, eval(ndim), evec, RLev_ENER(ndim),ZERO
      DOUBLE PRECISION A(ndim),B(ndim),GJ(ndim),SWAP1,SWAP2,SWAP3,SWAP4
!
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF, &
                     IMaxCount , MAX_STRING_LENGTH
      COMMON/HFS/ A,B,GJ
!
      DATA PlusMinus/'-', ' ', '+'/
      DATA JFraction/'  0 ', ' 1/2', '  1 ', ' 3/2', '  2 ', ' 5/2', &
                     '  3 ', ' 7/2', '  4 ', ' 9/2', '  5 ', '11/2', &
                     '  6 ', '13/2', '  7 ', '15/2', '  8 ', '17/2', &
                     '  9 ', '19/2', ' 10 ', &
                             '21/2', ' 11 ', '23/2', ' 12 ', '25/2', &
                     ' 13 ', '27/2', ' 14 ', '29/2', ' 15 ', '31/2', &
                     ' 16 ', '33/2', ' 17 ', '35/2', ' 18 ', '37/2', &
                     ' 19 ', '39/2', ' 20 ', '41/2', ' 21 ', '43/2', &
                     ' 22 '/
!
!
!     Opens the file  *.lsj.lbl
!

      WRITE(*,*)
      WRITE(*,*) 'RHFS_LSJ'
      WRITE(*,*) 'This program prints output from the rhfs program '
      WRITE(*,*) 'using LSJ lables. Output can be energy sorted'
      WRITE(*,*) 'Input files: name.(c)h, name.lsj.lbl'
      WRITE(*,*) 'Output file: name.(c)hlsj'
      WRITE(*,*)

      WRITE(*,*) 'Name of the state'
      READ(*,*) strFile 
      K = INDEX(strFile,' ')
      WRITE(*,*) 'Hfs data from a CI calc?'
      READ(*,*) ans
      IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
         hfs_file = strFile(1:K-1)//'.ch'
         output_file = strFile(1:K-1)//'.chlsj' 
      ELSE
         hfs_file = strFile(1:K-1)//'.h'
         output_file = strFile(1:K-1)//'.hlsj' 
      END IF
      OPEN (30, FILE = hfs_file, FORM = 'FORMATTED', &
         STATUS = 'OLD', IOSTAT = IOS)
      IF (IOS .NE. 0) THEN
         WRITE (0,*) 'Failed to open file ', TRIM(hfs_file)
         CLOSE (30)
         STOP
      END IF

      OPEN (80, FILE = output_file, FORM = 'FORMATTED', &
         STATUS = 'UNKNOWN')

      util_lbl_file = strFile(1:K-1)//'.lsj.lbl'
      OPEN (31, FILE = util_lbl_file, FORM = 'FORMATTED', &
         STATUS = 'OLD', IOSTAT = IOS)
      IF (IOS .NE. 0) THEN
         WRITE (0,*) 'Failed to open file ', TRIM(util_lbl_file)
         CLOSE (31)
         STOP
      END IF

      WRITE(*,*) 'Energy sorted output? '
      READ(*,*) ans
      IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
        nsort = 1
      ELSE
        nsort = 0
      END IF
       

      CALL READHFS
      CALL LDLBL

! Determine how many states there are with J = 0. These should not be
! printed since they have no hfs.

      IZERO = 0
      DO I = 1,IMaxCount
         IF (Lev_J(I).EQ.'   0') IZERO = IZERO + 1
      END DO       

      DUMMY0 = '                                '
      DUMMY = DUMMY0//DUMMY0
      N = MAX(1,MAX_STRING_LENGTH-4)
      WRITE(80,402) 'Energy',DUMMY(1:N),'State', &
         '    J   P          A(MHz)          B(MHz)         gJ'

      IF (nsort.eq.1) THEN
        DO J = IZERO + 2,IMaxCount
          SWAP1 = RLev_ENER(J)
          SWAP2 = A(J-IZERO)
          SWAP3 = B(J-IZERO)
          SWAP4 = GJ(J-IZERO)
          SWAPP = Lev_Par(J)
          SWAPJ = Lev_J(J)
          SWAPCSF = string_CSF(J)
          DO I = J-1,IZERO + 1,-1
            IF (RLev_ENER(I).LE.SWAP1) GOTO 10
            RLev_ENER(I+1) = RLev_ENER(I)
            A(I+1-IZERO) = A(I-IZERO)
            B(I+1-IZERO) = B(I-IZERO)
            GJ(I+1-IZERO) = GJ(I-IZERO)
            Lev_Par(I+1) = Lev_Par(I)
            Lev_J(I+1) = Lev_J(I)
            string_CSF(I+1) = string_CSF(I)
          END DO
          I = IZERO
   10     RLev_ENER(I+1) = SWAP1
          A(I+1-IZERO) = SWAP2
          B(I+1-IZERO) = SWAP3
          GJ(I+1-IZERO) = SWAP4
          Lev_Par(I+1) = SWAPP
          Lev_J(I+1) = SWAPJ
          string_CSF(I+1) = SWAPCSF
        END DO        
      END IF

      DO I = IZERO + 1,IMaxCount
!         WRITE(80,403) Lev_POS(I),Lev_J(I),Lev_Par(I),
!     &   A(I-IZERO),B(I-IZERO),GJ(I-IZERO),trim(string_CSF(I))
!         WRITE(80,404) Lev_POS(I),Lev_J(I),Lev_Par(I),
!     &   A(I-IZERO),B(I-IZERO),GJ(I-IZERO),trim(string_CSF(I)),
!     &   RLev_ENER(I)
         WRITE(80,405) RLev_ENER(I), &
         string_CSF(I)(1:MAX_STRING_LENGTH),Lev_J(I),Lev_Par(I), &
         A(I-IZERO),B(I-IZERO),GJ(I-IZERO)
      END DO

      STOP

  402 FORMAT (//9X,A,A,A,A)
!  402 FORMAT (//' Interaction constants:'//
!     :' Energy      J Parity ',1X,'A (MHz)',6X,'B (MHz)',8X,
!     : 'g_J            state'/)
  403 FORMAT (1X,1I3,5X,2A4,1P,2D13.3,D16.6,5X,A)
  404 FORMAT (1X,1I3,5X,2A4,1P,2D13.3,D16.6,5X,A,F14.7)
!  405 FORMAT (1X,F14.7,2X,A,2A4,1P,2D13.3,1P,D16.6)
  405 FORMAT (1X,F14.7,2X,A,2A4,3X,1P,D13.3,3X,1P,D13.3,1P,D16.6)


      END

!
!***********************************************************************
!                                                                      *
      SUBROUTINE READHFS
!                                                                      *
!     Open, check and load data from the hfs file                      *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by P. J\"onsson,                                         *
!     Malmo                                                 Aug 2011   *
!                                                                      *
!***********************************************************************
!
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER, PARAMETER :: ndim = 20000          ! max number of states

      CHARACTER*120 RECORD
      DOUBLE PRECISION A(ndim),B(ndim),GJ(ndim)
!
      COMMON/HFS/A,B,GJ 

! Position yourself at the correct place in the file

      DO I = 1,3
         READ (30,'(A)') RECORD
         WRITE(80,'(1X,A)') TRIM(RECORD)
      END DO
     

      DO
         READ (30,'(A)') RECORD
!PJ         WRITE(*,'(A)') RECORD
         IF (RECORD(30:32).EQ.'MHz') GOTO 10
      END DO
   10 CONTINUE
      READ (30,'(A)') RECORD

! Now begin to read the hfs data

      I = 1
      DO
         READ (30,'(17X,5D20.10)',IOSTAT=IOS) A(I),B(I),GJ(I)
         IF (IOS.NE.0) GOTO 20
!PJ         WRITE(*,'(1P,3D20.10)') A(I),B(I),GJ(I)
         I = I + 1
      END DO
   20 CONTINUE 

      RETURN
      END

!
!***********************************************************************
!                                                                      *
      SUBROUTINE LDLBL
!                                                                      *
!     Open, check and load data from the  .lsj.lbl   file of the       *
!     inital state.                                                    *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  May 2011   *
!                                                                      *
!***********************************************************************
!
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER, PARAMETER :: ndim = 20000          ! max number of states
!
      CHARACTER*1 Lev_par(ndim)
      CHARACTER*4 Lev_J(ndim)
      CHARACTER*15 RECORD
      CHARACTER*64 string_CSF(ndim)
!
      INTEGER IOS, ITEST, Lev_POS(ndim)
      INTEGER MAX_STRING_LENGTH
      REAL WEIGHTS
      DOUBLE PRECISION RLev_ENER(ndim)
!
      COMMON/JJ2LSJ/ Lev_POS,Lev_J,Lev_Par,RLev_ENER,string_CSF, &
                      IMaxCount, MAX_STRING_LENGTH
!
      MAX_STRING_LENGTH = 0 

      READ (31,'(1A15)',IOSTAT = IOS) RECORD
      ICount = 0
      IF (IOS .NE. 0) GO TO 1
      ICount = 1
      READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS) &
        Lev_Pos(ICount),Lev_J(ICount),Lev_Par(ICount), &
        RLev_ENER(ICount)
      IF (IOS .NE. 0) GO TO 1
!
      READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF(ICount)
      K = INDEX(string_CSF(ICount),' ')
      IF (K.GT.MAX_STRING_LENGTH) MAX_STRING_LENGTH = K
               
!
    2 READ (31,'(1X,I2)',IOSTAT = IOS) ITEST
      IF (IOS .NE. 0) GO TO 1
      IF (ITEST .EQ. 0) GO TO 2
      BACKSPACE 31
      ICount = ICount + 1
      READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS) &
        Lev_Pos(ICount),Lev_J(ICount),Lev_Par(ICount), &
        RLev_ENER(ICount)
      IF (IOS .NE. 0) GO TO 1
      READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF(ICount)
      K = INDEX(string_CSF(ICount),' ')
      IF (K.GT.MAX_STRING_LENGTH) MAX_STRING_LENGTH = K
      GO TO 2
    1 CONTINUE
      IMaxCount = ICount + IMaxCount
      RETURN
      END
