!                                                                      *
      SUBROUTINE OSCL(NAME)
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of  data for transitions between multiconfiguration   *
!   Dirac-Fock energy levels.                                          *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE blk_C
      USE biorb_C
      USE def_C, CCMPS=>CCMS
      USE default_C
      USE EIGV_C
      USE orb_C
      USE OSC_C
      USE PRNT_C
      USE SYMA_C
      USE TITL_C
      USE WAVE_C
      USE jj2lsjbio_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cpmix_I
      USE alcnsa_I
      USE alcnta_I
      USE connect_I
      USE mctout_I
      USE readmix_I
      USE mctin_I
      USE itrig_I
      USE bessj_I
      USE csfm_I
      USE printa_I
      USE printals_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*4, PARAMETER :: IAU = 'Hart'
      CHARACTER*4, PARAMETER :: IEV = ' eV '
      CHARACTER*4, PARAMETER :: ICM = 'Kays'
      CHARACTER*4, PARAMETER :: IHZ = ' Hz '
      CHARACTER*4, PARAMETER :: IANG = ' A  '
      INTEGER, PARAMETER :: NCA = 65536
      INTEGER, PARAMETER :: NFILE = 93
      INTEGER, PARAMETER :: NFILE1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KKA, JKP, NFILE2, IOS, my_NELEC, NCFTOTI, NWI, &
                 NVECSIZI, IOPAR, IBLKI, NCFTOTF, NWF,&
                 NVECSIZF, IBLKF, NVECPR, I, M,       &
                 IELEC, LEVII, LEVFF, ITKPO, ITEST, NLP, LINES, &
                 ILBL, INUM_II, INUM_FF, ICOUNT1, ICOUNT2
      REAL(DOUBLE) :: FACTOR, OMEGA, ARGU, ASFA, ASFB
      REAL(DOUBLE), DIMENSION(:), pointer :: et, et1, ipr, ipr1, next
      LOGICAL :: AVAIL, LSAME
      CHARACTER :: ANSW*1, IUNITS*4, G92MIX*6
!-----------------------------------------------
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list
!
! for the case that the inital file and final file are same
!
      LSAME = TRIM(NAME(1)) == TRIM(NAME(2))
      IF (LSAME) CALL CPMIX (NAME, INPCI)
!
! write header for the result file
!
      if (LSAME) then
         write(24,*) 'Transition in file:'
         write(24,*) 'f = ',trim(name(1))
      else
         write(24,*) 'Transition between files:'
         write(24,*) 'f1 = ',trim(name(1))
         write(24,*) 'f2 = ',trim(name(2))
      end if
!
      CALL ALCNSA (JJA, JJB, HB1, HB2, HC1, HC2, HM1, &
         HM2, LAB, NPTR, NSDIM, 1)
      CALL ALCNTA (ISLDR, ISLDR1, XSLDR, NTDIM, 1)
!
!   Make a connection between the orbitals of the
!   merged list and the initial and final state lists
!
      CALL CONNECT
!
!   Set up units for printing transition energy
!
      IF (LTC(1)) THEN
!
!   Print transition energies in Angstroms
!
         FACTOR = AUCM
         FACTOR = 1.0D08/FACTOR
         IUNITS = IANG
!
      ELSE IF (LTC(2)) THEN
!
!   Print energies in eV
!
         FACTOR = AUEV
         IUNITS = IEV
!
      ELSE IF (LTC(3)) THEN
!
!   Print transition energies in Hartree Atomic Units
!
         FACTOR = 1.0D00
         IUNITS = IAU
!
      ELSE IF (LTC(4)) THEN
!
!   Print transition energies in Hz
!
         FACTOR = AUCM
         FACTOR = FACTOR*CCMPS
         IUNITS = IHZ
!
      ELSE IF (LTC(5)) THEN
!
!   Print transition energies in Kaysers
!
         FACTOR = AUCM
         IUNITS = ICM
!
      ENDIF
!
!   Select type of transition
!
!   KK  =  0  for electric multipole.
!       =  1  for magnetic multipole.
!
!   CF:  IOPAR  =  (-1)**N      Electric N-pole
!               =  (-1)**(N+1)  Magnetic N-pole.
!                               N > 0
!
      KKA = 1
      DO JKP = 1, NKP
         NFILE2 = NFILE1 + JKP
!
! read the head of the file of mixing coef. for initial
!
         IF (LSAME) NAME(1) = TRIM(NAME(2))//'_CP'
         IF (INPCI == 0) THEN
            OPEN(UNIT=68, FILE=TRIM(NAME(1))//'.cbm', FORM='UNFORMATTED', &
              STATUS='OLD')
         ELSE
            OPEN(UNIT=68, FILE=TRIM(NAME(1))//'.bm', FORM='UNFORMATTED', &
               STATUS='OLD')
         ENDIF
         IF (LSAME) NAME(1) = NAME(2)
         READ (68, IOSTAT=IOS) G92MIX
         IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
            WRITE (*, *) 'Not a GRASP mixing file'
            STOP
         ENDIF
         READ (68) NELEC, NCFTOTI, NWI, NVECTOTI, NVECSIZI, NBLOCKI
         WRITE (*, *) '   nelec  = ', my_NELEC
         WRITE (*, *) '   ncftoti = ', NCFTOTI
         WRITE (*, *) '   nwi     = ', NWI
         WRITE (*, *) '   nblocki = ', NBLOCKI
         WRITE (*, *)
!   lbl
         ICOUNT1 = 0
         CALL LDLBL1 (NAME(1))
! If not available generate angular coefficients for all pares of blocks
         CALL MCTOUT (IOPAR, JKP, NAME)
         DO IBLKI = 1, NBLOCKI
            CALL READMIX (NAME, INPCI, 1)
!   lbl
            ICOUNT1 = NVECII + ICOUNT1
!
! read the head of the file of mixing coef. for final
!
            IF (INPCI == 0) THEN
               OPEN(UNIT=78,FILE=TRIM(NAME(2))//'.cbm',              &
                  FORM='UNFORMATTED',STATUS='OLD')
            ELSE
               OPEN(UNIT=78, FILE=TRIM(NAME(2))//'.bm',              &
                  FORM='UNFORMATTED',STATUS='OLD')
            ENDIF
            READ (78, IOSTAT=IOS) G92MIX
            IF (IOS/=0 .OR. G92MIX/='G92MIX') THEN
               WRITE (*, *) 'Not a GRASP mixing file'
               STOP
            ENDIF
            READ (78) NELEC, NCFTOTF, NWF, NVECTOTF, NVECSIZF, NBLOCKF
            WRITE (*, *) '   nelec  = ', my_NELEC
            WRITE (*, *) '   ncftotf = ', NCFTOTF
            WRITE (*, *) '   nwf     = ', NWF
            WRITE (*, *) '   nblockf = ', NBLOCKF
            WRITE (*, *)
!GG  lbl
            IF(IBLKI .EQ. 1) THEN
               CALL LDLBL2 (NAME(2))
               IF(IOPEN_STATUS1.EQ.0 .AND. IOPEN_STATUS2 .EQ.0) THEN
                  ILBL = 1
                  CALL STRSUM(NAME,INPCI,ILBL)
                  write(32,*) 'Transition between files:'
                  write(32,*) trim(name(1))
                  write(32,*) trim(name(2))
               END IF
            END IF
            ICOUNT2 = 0
!GG   lbl end
            DO IBLKF = 1, NBLOCKF
               CALL READMIX (NAME, INPCI, 2)
!   lbl
               ICOUNT2 = NVECFF + ICOUNT2
!
!   Allocate storage
!
               CALL ALLOC (TOTB, NVECFF, 'TOTB', 'OSCL')
               CALL ALLOC (TOTC, NVECFF, 'TOTC', 'OSCL')
!
               NVECPR = NVECII*NVECFF
               CALL ALLOC (ET, NVECPR, 'ET', 'OSCL' )
               CALL ALLOC (ET1, NVECPR, 'ET1','OSCL' )
               CALL ALLOC (IPR, NVECPR, 'IPR', 'OSCL' )
               CALL ALLOC (IPR1, NVECPR,'IPR1', 'OSCL' )
               CALL ALLOC (NEXT, NVECPR, 'NEXT', 'OSCL')
!
!   Initialization for total decay rate
!
               TOTC(:NVECFF) = 0.0D00
               TOTB(:NVECFF) = 0.0D00
!
               CALL MCTIN (IOPAR, JKP, NAME)
!
               IF (LK > 0) THEN
                  IELEC = (-1)**LK
                  IF (IELEC == IOPAR) THEN
                     KK = 0
                     KKA = 0
                  ELSE
                     KK = 1
                     IELEC = -IELEC
                  ENDIF
!
!   Set up list of levels for calculation of oscillator strengths
!   sort list into increasing order of energy if option 6 set
!
                  IF (IBLKI==1 .AND. IBLKF==1) THEN
                     IF (KK == 0) THEN
                        WRITE (24, 308) LK
                     ELSE
                        WRITE (24, 309) LK
                     ENDIF
                     WRITE (24, 310)
                     IF (LTC(1)) THEN
                        WRITE (24, 311)
                        IF (.NOT.LTC(7)) THEN
                           WRITE (24, 312)
                        ELSE
                           WRITE (24, 313)
                        ENDIF
                     ELSE
                        WRITE (24, 314)
                        IF (.NOT.LTC(7)) THEN
                           WRITE (24, 315) IUNITS
                        ELSE
                           WRITE (24, 316) IUNITS
                        ENDIF
                     ENDIF
                  ENDIF
!
                  DO LEVII = 1, NVECII
!   lbl
                     INUM_II = ICOUNT1 - NVECII + LEVII
                     DO LEVFF = 1, NVECFF
!   lbl
                        INUM_FF = ICOUNT2 - NVECFF + LEVFF
!
!   Check for consistent parity and J
!
                        ITKPO = LK + LK + 1
                        IF (ITRIG(IATJPOII(LEVII),IATJPOFF(LEVFF),ITKPO) == 0) &
                           CYCLE
                        ITEST = IASPARII(LEVII)*IASPARFF(LEVFF)*IELEC
                        IF (ITEST < 0) CYCLE
!
!   Calculate and print transition probability data
!
                        NLP = 70 - 8
                        LINES = NLP
!
                        IF (LINES >= NLP) LINES = 0
!
                        M = LEVFF + NVECFF*(LEVII - 1)
!             M = LEVFF+NVECII*(LEVII-1)
                        ET(M) = EVALFF(LEVFF) + EAVFF - EVALII(LEVII) - EAVII
                        IF (LSAME .AND. ET(M)<=0.0) CYCLE
                        OMEGA = -ET(M)
                        ARGU = OMEGA/C
                        CALL BESSJ (ARGU)
!
!  Calculate oscillator strength between the ASFs
!
                        CALL CSFM (ASFA, ASFB, LEVII, LEVFF)
                        CALL PRINTA (ASFA, ASFB, LEVII, LEVFF, OMEGA, FACTOR, &
                           LINES,LSAME)
               IF(IOPEN_STATUS1.EQ.0 .AND. IOPEN_STATUS2 .EQ.0) THEN
                 CALL PRINTALS (INUM_II, INUM_FF,ASFA,ASFB,LEVII,LEVFF,&
                           OMEGA,FACTOR)
               END IF

!              WRITE (24,317)
                     END DO
                  END DO
               ENDIF
!
!   Deallocate storage; this is local to OSCL
!
               CALL DALLOC (TOTB, 'TOTB', 'OSCL')
               CALL DALLOC (TOTC, 'TOTC', 'OSCL')
!
               CALL DALLOC (ET, 'ET', 'OSCL')
               CALL DALLOC (ET1, 'ET!', 'OSCL')
               CALL DALLOC (IPR, 'IPR', 'OSCL')
               CALL DALLOC (IPR1, 'IPR1', 'OSCL')
               CALL DALLOC (NEXT, 'NEXT', 'OSCL')
!
!   Deallocate storage; this is allocated in READMIX
!
               CALL DALLOC (EVALFF, 'EVALFF', 'OSCL')
               CALL DALLOC (EVECFF, 'EVECFF', 'OSCL')
               CALL DALLOC (IVECFF, 'IVECFF', 'OSCL')
               CALL DALLOC (IATJPOFF, 'IATJPOFF', 'OSCL')
               CALL DALLOC (IASPARFF, 'IASPARFF', 'OSCL')

            END DO
            CLOSE(78)
!
!   Deallocate storage; this is allocated in READMIX
!
            CALL DALLOC (EVALII, 'EVALII', 'OSCL')
            CALL DALLOC (EVECII, 'EVECII', 'OSCL')
            CALL DALLOC (IVECII, 'IVECII', 'OSCL')
            CALL DALLOC (IATJPOII, 'IATJPOII', 'OSCL')
            CALL DALLOC (IASPARII, 'IASPARII', 'OSCL')

         END DO
         CLOSE(68)
         CLOSE(NFILE2)
      END DO
      CALL DALLOC (KP, 'KP', 'OSCL')
!
! close and delete duplicated mixing file
      IF (LSAME) THEN
         IF (INPCI == 0) THEN
!GG            OPEN(68, FILE=TRIM(NAME(2))//'_CP.cbm', POSITION='asis')
            OPEN(68, FILE=TRIM(NAME(2))//'_CP.cbm')
         ELSE
!GG            OPEN(68, FILE=TRIM(NAME(2))//'_CP.bm', POSITION='asis')
            OPEN(68, FILE=TRIM(NAME(2))//'_CP.bm')
         ENDIF
         CLOSE(68, STATUS='delete')
      ENDIF
!
      CALL ALCNSA (jja, jjb, hb1, hb2, hc1, hc2, hm1, hm2, lab, nptr, nsdim, 3)
      CALL ALCNTA (isldr, isldr1, xsldr, ntdim, 3)
!
!
!   Close all files
!
      CLOSE(24)
!
      RETURN
!
  302 FORMAT(/,' ***** Warning *****')
  303 FORMAT(/,/,' ***** Error in OSCL *****')
  307 FORMAT(/,' Dynamic allocation computed incorrectly: Bug.')
  308 FORMAT(/,/,' Electric 2**(',I2,')-pole transitions')
  309 FORMAT(/,/,' Magnetic 2**(',I2,')-pole transitions')
  310 FORMAT(1X,33('='))
  311 FORMAT(/,'   Upper state        Lower state  ',8X,'Gauge',8X,'Wavelength'&
         ,13X,'Einstein coefficients',13X,'Oscillator')
  312 FORMAT(81X,'-1',15X,'3 -2 -1',/,' Level  J Parity',4X,'Level  J ',&
         'Parity',21X,'(Angstroms)',10X,'A (s  )',9X,'gB (m s  J  )',7X,&
         'strength gf'/)
  313 FORMAT(' Level  J Parity',4X,'Level  J Parity',21X,'(Angstroms)',10X,&
         'A (au)',13X,'gB (au)',10X,'strength gf'/)
  314 FORMAT(/,' Upper       Lower ')
  315 FORMAT(' Lev  J P',3X,'Lev  J P',7X,'E (',A4,')',9X,'A (s-1)',10X,'gf',&
         12X,'S', 13X,'M')
  316 FORMAT(' Level  J Parity',4X,'Level  J Parity',23X,'(',A4,')',13X,&
         'A (au)',13X,'gB (au)',10X,'strength gf'/)
  317 FORMAT(/,1X,124('+'))
  318 FORMAT(/,/,' Radiative lifetimes '/,' ======================='/,/,&
         ' Level      Lifetime  s (-1)')
  319 FORMAT(1X,I4,6X,'Coulomb: ',1P,1D20.7)
  320 FORMAT(10X,'Babushkin:',1P,1D20.7,/)
  321 FORMAT(1X,I4,5X,'Magnetic: ',1P,1D20.7,/)
      RETURN
!
      END SUBROUTINE OSCL
