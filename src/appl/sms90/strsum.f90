!***********************************************************************
!                                                                      *
      SUBROUTINE STRSUM 
!                                                                      *
!   Generates the first part of  sms92.sum  (on stream 24).            *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, CONVRT, WGHTD5.                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE decide_C
      USE def_C
      USE eigv_C
      USE iccu_C
      USE grid_C
      USE npar_C
      USE npot_C
      USE nsmdat_C
      USE orb_C
      USE prnt_C
      USE syma_C
      USe wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I 
      USE engout_I 
      USE wghtd5_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: LENTH, I 
      CHARACTER :: RECORD*256, CDATA*26 
!-----------------------------------------------
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
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
         CALL CONVRT (ICCUT(1), RECORD, LENTH) 
         WRITE (24, *) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'//&
            ' the zero-order space;' 
      ENDIF 
!
!   Write out the nuclear parameters
!
      WRITE (24, *) 
      WRITE (24, 300) Z 
      IF (NPARM == 2) THEN 
         WRITE (24, *) 'Fermi nucleus:' 
         WRITE (24, 301) PARM(1), PARM(2) 
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
      WRITE (24, 305) C 
!
!   Write out the parameters of the radial grid
!
      WRITE (24, *) 
      IF (HP == 0.0D00) THEN 
         WRITE (24, 306) RNT, H, N 
      ELSE 
         WRITE (24, 307) RNT, H, HP, N 
      ENDIF 
      WRITE (24, 308) R(1), R(2), R(N) 
!
!   Write out the orbital properties
!
      WRITE (24, *) 
      WRITE (24, *) 'Subshell radial wavefunction summary:' 
      WRITE (24, *) 
      WRITE (24, 309) 
      WRITE (24, *) 
      DO I = 1, NW 
         WRITE (24, 310) NP(I), NH(I), E(I), PZ(I), GAMA(I), PF(2,I), QF(2,I), &
            MF(I) 
      END DO 
!
!   Write the list of eigenpair indices
!
      WRITE (24, *) 
      CALL ENGOUT (EAV, EVAL, IATJPO, IASPAR, IVEC, NVEC, 3) 
      CALL WGHTD5 
!
      RETURN  
!
  300 FORMAT('The atomic number is ',1F14.10,';') 
  301 FORMAT('  c =',1P,1D19.12,' Bohr radii,'/,'  a =',1D19.12,' Bohr radii;') 
  305 FORMAT('Speed of light = ',1P,D19.12,' atomic units.') 
  306 FORMAT('Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',' I = 1, ..., N;'/,/,&
         ' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D19.12,' Bohr radii;'/&
         ,' N    = ',1I4,';') 
  307 FORMAT('Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
         ' I = 1, ..., N;'/,/,' RNT  = ',1P,D19.12,' Bohr radii;'/,' H    = ',D&
         19.12,' Bohr radii;'/,' HP   = ',D19.12,' Bohr radii;'/,' N    = ',1I4&
         ,';') 
  308 FORMAT(' R(1) = ',1P,1D19.12,' Bohr radii;'/,' R(2) = ',1D19.12,&
         ' Bohr radii;'/,' R(N) = ',1D19.12,' Bohr radii.') 
  309 FORMAT(' Subshell',11X,'e',20X,'p0',18X,'gamma',19X,'P(2)',18X,'Q(2)',10X&
         ,'MTP') 
  310 FORMAT(3X,1I2,1A2,1X,1P,5(3X,1D19.12),3X,1I3) 
      RETURN  
!
      END SUBROUTINE STRSUM 
