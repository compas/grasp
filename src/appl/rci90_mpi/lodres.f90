!***********************************************************************
!                                                                      *
      SUBROUTINE LODRES 
!                                                                      *
!   Loads the data from the  .res  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!   This subroutine is called by every node and the only differences   *
!   from the non-mpi version are:                                      *
!     1) STOP and stopmpi (not necessary)                              *
!     2) LSE needs coomunication - can be moved out of this routine.   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Block version by Xinghong He          Last revision:  1 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNP
      USE memory_man
      USE decide_C
      USE def_C, ONLY: nelecr, nelec, z, emn, c, accy
      USE grid_C
      USE npar_C
      USE npot_C, ONLY: zz, nnuc
      USE orb_C
      USE wave_C
      USE wfac_C
      USE where_C
      USE hblock_C
      USE iccu_C
      USE iounit_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE setqic_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NCFRES, NWRES, NBLOCKRES, I, NP10, J 
      LOGICAL :: YES 
!-----------------------------------------------
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))
!
!-----------------------------------------------------------------------
!
!   Read the basic parameters of the electron cloud; check these
!   against those deduced from the  .csl file
!
      READ (IMCDF) NELECR, NCFRES, NWRES, NBLOCKRES 
 
      IF (NELECR/=NELEC .OR. NCFRES/=NCF .OR. NWRES/=NW .OR. NBLOCKRES/=NBLOCK&
         ) CALL STOPMPI ('lodres: NELEC/NCF/NW does not match',myid)
!
!   Read the nuclear parameters
!
      READ (IMCDF) Z, EMN 
      READ (IMCDF) NPARM, (PARM(I),I=1,NPARM) 
      READ (IMCDF) N, (ZZ(I),I=1,N), NNUC 
 
      IF (N > NNNP) CALL STOPMPI ('lodres: N greater than NNNP',myid)
!
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
      READ (IMCDF) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, WFACT, LVP, &
         LNMS, LSMS 
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
      NP10 = N + 10 
      READ (IMCDF) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), (RPOR(I),I=1,&
         NP10) 
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6 
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC 
!
!   Allocate storage for the radial wavefunction arrays
!
      CALL ALLOC (PF, NNNP,NW, 'PF', 'LODMIX') 
      CALL ALLOC (QF, NNNP,NW, 'QF', 'LODMIX') 
!
!   Read the orbital wavefunctions and the associated arrays
!
      DO J = 1, NW 
         READ (IMCDF) E(J), GAMA(J), PZ(J), MF(J) 
         READ (IMCDF) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J)) 
      END DO 
!
!   Determine if the self-energy contribution is to be estimated
!
      IF (MYID == 0) THEN
         WRITE (ISTDE, *) 'Estimate contributions from the self-energy?'
         LSE = GETYN() 
         IF (LSE) THEN
           WRITE(734,'(a)') 'y            ! Estimate contributions from the self-energy?'
         ELSE
           WRITE(734,'(a)') 'n            ! Estimate contributions from the self-energy?'
         END IF
      ENDIF
      CALL MPI_BCAST (LSE, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR) 
      RETURN  
      END SUBROUTINE LODRES 
