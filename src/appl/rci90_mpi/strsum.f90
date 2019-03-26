!***********************************************************************
!                                                                      *
      SUBROUTINE STRSUM
!                                                                      *
!   Generates the first part of  grasp92.sum  (on stream 24).          *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, convrt2.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 09 Dec 1992   *
!   Modified by Xinghong He               Last revision: 22 Dec 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      Use decide_C
      USE def_C
      USE grid_C
      USE npar_C
      USE npot_C
      USE orb_C
      USE prnt_C
      USE wave_C
      USE wfac_C
      USE iounit_C
      USE hblock_C
      USE iccu_C,          ONLY: ICCUTBLK
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt2_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: LENTH, NB, ICCUT, I, NTMP
      CHARACTER :: RECORD*256, CDATA*26, CLEVEL*2
!-----------------------------------------------
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))

!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!   Write out the basic dimensions of the electron cloud
!
      WRITE (24, *)
      CALL CONVRT2 (NELEC, RECORD, LENTH, 'strsum.nelec')
      WRITE (24, *) 'There are '//RECORD(1:LENTH)//' electrons in the cloud'
      CALL CONVRT2 (NCF, RECORD, LENTH, 'strsum.ncf')
      WRITE (24, *) ' in '//RECORD(1:LENTH)//' relativistic CSFs'
      CALL CONVRT2 (NW, RECORD, LENTH, 'strsum.nw')
      WRITE (24, *) ' based on '//RECORD(1:LENTH)//' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      IF (LFORDR) THEN
         WRITE (24, *)
         DO NB = 1, NBLOCK
            ICCUT = ICCUTBLK(NB)
            CALL CONVRT2 (ICCUT, RECORD, LENTH, 'strsum.icccut')
            WRITE (24, *) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'//&
               ' the zero-order space;  nb = ', NB, ' ncf = ', NCFBLK(NB)
         END DO
      ENDIF
!
!   Write out the nuclear parameters
!
      WRITE (24, *)
      WRITE (24, 300) Z
      IF (EMN == 0.0D00) THEN
         WRITE (24, *) ' the nucleus is stationary;'
      ELSE
         WRITE (24, 301) EMN
      ENDIF
      IF (NPARM == 2) THEN
         WRITE (24, *) ' Fermi nucleus:'
         WRITE (24, 302) PARM(1), PARM(2)
         CALL CONVRT2 (NNUC, RECORD, LENTH, 'strsum.nnuc')
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
      WRITE (24, *)
      IF (LTRANS .OR. LVP .OR. LNMS .OR. LSMS) THEN
         WRITE (24, *) 'To H (Dirac Coulomb) is added'
         IF (LTRANS) WRITE (24, 304) WFACT
         IF (LVP) WRITE (24, *) ' H (Vacuum Polarisation);'
         IF (LNMS) WRITE (24, *) ' H (Normal Mass Shift);'
         IF (LSMS) WRITE (24, *) ' H (Specific Mass Shift);'
         WRITE (24, *) ' the total will be diagonalised.'
      ELSE
         WRITE (24, *) 'H (Dirac Coulomb) will be diagonalised by itself.'
      ENDIF
!
      IF (LSE) THEN
         WRITE (24, *) &
            'Diagonal contributions from H (Self Energy) will be estimated'
         WRITE (24, *) ' from a screened hydrogenic approximation.'
      ENDIF
!
!   Write out the parameters of the radial grid
!
      WRITE (24, *)
      IF (HP == 0.0D00) THEN
         WRITE (24, 305) RNT, H, N
      ELSE
         WRITE (24, 306) RNT, H, HP, N
      ENDIF
      WRITE (24, 307) R(1), R(2), R(N)
!
!   Write out the orbital properties
!
      WRITE (24, *)
      WRITE (24, *) 'Subshell radial wavefunction summary:'
      WRITE (24, *)
      WRITE (24, 308)
      WRITE (24, *)
      DO I = 1, NW
         WRITE (24, 309) NP(I), NH(I), E(I), PZ(I), GAMA(I), PF(2,I), QF(2,I), &
            MF(I)
      END DO
!
!   Write the list of eigenpair indices
!
      WRITE (24, *)
!
!   Find total number of eigenstates and print corresponding info
!
      ntmp = 0
      DO i = 1, nblock
         ntmp = ntmp + nevblk(i)
      ENDDO

      WRITE (24,*) ntmp, ' levels will be computed'

      RETURN
!
  300 FORMAT('The atomic number is ',1F14.10,';')
  301 FORMAT(' the mass of the nucleus is ',1P,D19.12,' electron masses;')
  302 FORMAT('  c =',1P,1D19.12,' Bohr radii,'/,'  a =',1D19.12,' Bohr radii;')
  303 FORMAT('Speed of light = ',1P,D19.12,' atomic units.')
  304 FORMAT('  H (Transverse) --- factor multiplying the',&
         ' photon frequency: ',1P,D15.8,';')
  305 FORMAT('Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',' I = 1, ..., N;'/,/,&
         ' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D19.12,' Bohr radii;'/&
         ,' N    = ',1I4,';')
  306 FORMAT('Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
         ' I = 1, ..., N;'/,/,' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D&
         19.12,' Bohr radii;'/,' HP   = ',D19.12,' Bohr radii;'/,' N    = ',1I4&
         ,';')
  307 FORMAT(' R(1) = ',1P,1D19.12,' Bohr radii;'/,' R(2) = ',1D19.12,&
         ' Bohr radii;'/,' R(N) = ',1D19.12,' Bohr radii.')
  308 FORMAT('Subshell',6X,'e',13X,'p0',5X,'gamma',5X,'P(2)',7X,'Q(2)',4X,'MTP'&
         )
  309 FORMAT(1X,I2,A2,1X,1P,D17.10,1P,D11.3,0P,F6.2,1P,2(D11.3),I5)
      RETURN
!
      END SUBROUTINE STRSUM
