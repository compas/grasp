!***********************************************************************
!                                                                      *
      SUBROUTINE STRSUM(NAME, INPCI, ILBL)
!                                                                      *
!   Generates the first part of  oscl92.sum  (on stream 24).           *
!                                                                      *
!   Call(s) to: [LIB92]: CALEN, CONVRT, WGHTD5.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE biorb_C
      USE decide_C
      USE def_C
      USE foparm_C
      USE grid_C
      USE eigv_C
      USE npar_C
      USE orb_C
      USE prnt_C
      USE syma_C
      USE wave_C
      USE jj2lsjbio_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: INPCI
      INTEGER :: ILBL
      CHARACTER, INTENT(IN) :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J
      CHARACTER :: RECORD*15, CTIME*8, CDATE*8
!-----------------------------------------------


      I = INDEX(NAME(1),' ')
      J = INDEX(NAME(2),' ')
      IF(ILBL == 0) THEN
         IF (INPCI == 0) THEN
            OPEN(UNIT=24,FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.ct',FORM=&
            'FORMATTED', STATUS='UNKNOWN',POSITION='asis')
         ELSE
            OPEN(UNIT=24,FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.t',FORM=&
            'FORMATTED',STATUS='UNKNOWN',POSITION='asis')
         ENDIF
      ELSE IF(ILBL == 1) THEN
         IF(IOPEN_STATUS1.EQ.0 .AND. IOPEN_STATUS2 .EQ.0) THEN
            IF (INPCI == 0) THEN
               OPEN(UNIT=32,                                           &
                   FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.ct.lsj',&
                   FORM='FORMATTED',STATUS='UNKNOWN')
            ELSE
               OPEN(UNIT=32,                                           &
                   FILE=NAME(1)(1:I-1)//'.'//NAME(2)(1:J-1)//'.t.lsj', &
                   FORM='FORMATTED',STATUS='UNKNOWN')
            ENDIF
         END IF
      END IF
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!ww      CALL CALEN (CTIME,CDATE)
!ww      WRITE (24,*) 'OSCL92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the initial state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECII,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the initial state cloud'
!ww      CALL CONVRT (NCFII,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWII,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   Write out the basic dimensions of the final state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECFF,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the final state cloud'
!ww      CALL CONVRT (NCFFF,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWFF,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
!ww      IF (LFORDR) THEN
!ww         WRITE (24,*)
!ww         CALL CONVRT (ICCUT,RECORD,LENTH)
!ww         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
!ww     :              //' the zero-order space;'
!ww      ENDIF
!
!   Write out the nuclear parameters
!
!ww      WRITE (24,*)
!ww      WRITE (24,300) Z
!ww      IF (NPARM .EQ. 2) THEN
!ww         WRITE (24,*) 'Fermi nucleus:'
!ww         WRITE (24,301) PARM(1),PARM(2)
!ww      ELSE
!ww         WRITE (24,*) ' point nucleus.'
!ww      ENDIF
!
!   Write out the physical effects specifications
!
!ww      WRITE (24,*)
!ww      WRITE (24,305) C
!
!   Write out the parameters of the radial grid
!
!ww      WRITE (24,*)
!ww      IF (HP .EQ. 0.0D 00) THEN
!ww         WRITE (24,306) RNT,H,N
!ww      ELSE
!ww         WRITE (24,307) RNT,H,HP,N
!ww      ENDIF
!ww      WRITE (24,308) R(1),R(2),R(N)
!
!   Write out the orbital properties
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Initial state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 1 I = 1,NWII
!ww         WRITE (24,310) NPII(I),NHII(I),EII(I),PZII(I),
!ww     :                  GAMAII(I),PFII(2,I),QFII(2,I),MFII(I)
!ww    1 CONTINUE
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Final state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 2 I = 1,NWFF
!ww         WRITE (24,310) NPFF(I),NHFF(I),EFF(I),PZFF(I),
!ww     :                  GAMAFF(I),PFFF(2,I),QFFF(2,I),MFFF(I)
!ww    2 CONTINUE
!
!   Write the list of eigenpair indices for the initial state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVII,EVALII,IATJPOII,IASPARII,IVECII,NVECII,3,1)
!
!   Write the list of eigenpair indices for the final state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVFF,EVALFF,IATJPOFF,IASPARFF,IVECFF,NVECFF,3,2)
!
      RETURN
!
  300 FORMAT ('The atomic number is ',1F14.10,';')
  301 FORMAT ('  c =',1P,1D19.12,' Bohr radii,'            &
             /'  a =',   1D19.12,' Bohr radii;')
  305 FORMAT ('Speed of light = ',1PD19.12,' atomic units.')
  306 FORMAT ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',&
               ' I = 1, ..., N;'                           &
             //' RNT  = ',1P,D19.12,' Bohr radii;'         &
              /' H    = ',   D19.12,' Bohr radii;'         &
              /' N    = ',1I4,';')
  307 FORMAT ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
               ' I = 1, ..., N;'                           &
             //' RNT  = ',1P,D19.12,' Bohr radii;'         &
              /' H    = ',   D19.12,' Bohr radii;'         &
              /' HP   = ',   D19.12,' Bohr radii;'         &
              /' N    = ',1I4,';')
  308 FORMAT ( ' R(1) = ',1P,1D19.12,' Bohr radii;'        &
              /' R(2) = ',   1D19.12,' Bohr radii;'        &
              /' R(N) = ',   1D19.12,' Bohr radii.')
  309 FORMAT (' Subshell',11X,'e',20X,'p0',18X,            &
              'gamma',19X,'P(2)',18X,'Q(2)',10X,'MTP')
  310 FORMAT (3X,1I2,1A2,1X,1P,5(3X,1D19.12),3X,1I3)

      RETURN
!
      END SUBROUTINE STRSUM
