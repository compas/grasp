!***********************************************************************
      SUBROUTINE STRSUM

!   Generates the first part of  rscf92.sum  (on stream 24).
!
!   Call(s) to: [LIB92] CALEN, CONVRT.
!
!   Written by Farid A. Parpia            Last revision: 26 Sep 1993
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE foparm_C
      USE grid_C
      USE mcpb_C
      USE npar_C
      USE npot_C
      USE orb_C
      USE wave_C
      USE wfac_C
!      USE iccu_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENTH, IEND, I, IBEG
      CHARACTER :: RECORD*256, CDATA*26, CTIME*8, CDATE*8, CLEVEL*2
!-----------------------------------------------
!
!     POINTER (PWEIGH,WEIGHT(1))
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRPF,PF(NNNP,1))
!     POINTER (PNTRQF,QF(NNNP,1))
!
!   Both the nuclear charge and the number of electrons are
!   known at this point; load IONCTY with the ionicity
!
      IONCTY = NINT(Z) - NELEC
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!      CALL CALEN (CTIME, CDATE)
!      WRITE (24,*) 'RSCF92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the electron cloud
!
      WRITE (24, *)
      CALL CONVRT (NELEC, RECORD, LENTH)
      WRITE (24, *) 'There are '//RECORD(1:LENTH)//' electrons in the cloud'
      CALL CONVRT (NCF, RECORD, LENTH)
      WRITE (24, *) ' in '//RECORD(1:LENTH)//' relativistic CSFs'
      CALL CONVRT (NW, RECORD, LENTH)
      WRITE (24, *) ' based on '//RECORD(1:LENTH)//' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      IF (LFORDR) THEN
         WRITE (24, *)
         CALL CONVRT (ICCUT, RECORD, LENTH)
         WRITE (24, *) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'//&
            ' the zero-order space;'
      ENDIF
!
!   Write out the nuclear parameters
!
      WRITE (24, *)
      WRITE (24, 300) Z
      IF (EMN == 0.D0) THEN
         WRITE (24, *) ' the nucleus is stationary;'
      ELSE
         WRITE (24, 301) EMN
      ENDIF
      IF (NPARM == 2) THEN
         WRITE (24, *) ' Fermi nucleus:'
         WRITE (24, 302) PARM(1), PARM(2)
         CALL CONVRT (NNUC, RECORD, LENTH)
         WRITE (24, *) ' there are '//RECORD(1:LENTH)//&
            ' tabulation points in the nucleus.'
      ELSE
         WRITE (24, *) ' point nucleus.'
      ENDIF
!
!   Write out the physical effects specifications
!
      WRITE (24, *)
      WRITE (24, 303) C
!
!   Write out the parameters of the radial grid
!
      WRITE (24, *)
      IF (HP == 0.D0) THEN
         WRITE (24, 305) RNT, H, N
      ELSE
         WRITE (24, 306) RNT, H, HP, N
      ENDIF
      WRITE (24, 307) R(1), R(2), R(N)
      WRITE (24, *)
!
!  (E)AL calculation, returns here
!
      IF (NCMIN == 0) THEN
         WRITE (24, *) '(E)AL calculation.'
         RETURN
      ENDIF
!
!  Info exclusively for EOL calculations
!
      IF (NCMIN == 1) THEN
         WRITE (24, *) 'OL calculation.'
         CALL CONVRT (ICCMIN(1), RECORD, LENTH)
         WRITE (24, *) 'Level '//RECORD(1:LENTH)//' will be optimised.'
      ELSE
         WRITE (24, *) 'EOL calculation.'
         CALL CONVRT (NCMIN, RECORD, LENTH)
         WRITE (24, *) RECORD(1:LENTH)//' levels will be optimised;'
         RECORD(1:20) = ' their indices are: '
         IEND = 20
         DO I = 1, NCMIN
            IBEG = IEND + 1
            CALL CONVRT (ICCMIN(I), CLEVEL, LENTH)
            IF (I /= NCMIN) THEN
               IEND = IBEG + LENTH + 1
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//', '
            ELSE
               IEND = IBEG + LENTH
               RECORD(IBEG:IEND) = CLEVEL(1:LENTH)//'.'
            ENDIF
            IF (IEND < 120) CYCLE
            WRITE (24, *) RECORD(1:IEND)
            RECORD(1:2) = '  '
            IEND = 2
         END DO
         IF (IEND /= 2) WRITE (24, *) RECORD(1:IEND)
         IF (WEIGHT(1) == (-1.D0)) THEN
            WRITE (24, *) 'Each is assigned its statistical weight;'
         ELSE IF (WEIGHT(1) == (-2.D0)) THEN
            WRITE (24, *) 'All levels are weighted equally;'
         ELSE
            WRITE (24, *) ' weighted as follows:'
            WRITE (24, *) (WEIGHT(I),I=1,NCMIN)
         ENDIF
      ENDIF

  300 FORMAT('The atomic number is ',1F14.10,';')
  301 FORMAT(' the mass of the nucleus is ',1P,D19.12,' electron masses;')
  302 FORMAT('  c =',1P,1D19.12,' Bohr radii,'/,'  a =',1D19.12,' Bohr radii;')
  303 FORMAT('Speed of light = ',3P,D19.12,' atomic units.')
  305 FORMAT('Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',' I = 1, ..., N;'/,/,&
         ' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D19.12,' Bohr radii;'/&
         ,' N    = ',1I4,';')
  306 FORMAT('Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
         ' I = 1, ..., N;'/,/,' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D&
         19.12,' Bohr radii;'/,' HP   = ',D19.12,' Bohr radii;'/,' N    = ',1I4&
         ,';')
  307 FORMAT(' R(1) = ',1P,1D19.12,' Bohr radii;'/,' R(2) = ',1D19.12,&
         ' Bohr radii;'/,' R(N) = ',1D19.12,' Bohr radii.')

      RETURN
      END SUBROUTINE STRSUM
