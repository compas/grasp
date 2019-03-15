!***********************************************************************
      SUBROUTINE IMPROVmpi (EOL, J, LSORT, DAMPMX)
!   The difference from the serial version is that it calls MPI        *
!   version subroutines (setlagmpi, cofpotmpi, matrixmpi, newcompi).   *
!                                                                      *
!   Improve the orbital J.                                             *
!                                                                      *
!   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, matrixmpi,    *
!                         newcompi, ORTHOR, ROTATE, SETCOF, setlagmpi, *
!                         SOLVE, XPOT, YPOT.                           *
!               [LIB92]: ORTHSC, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   Modified by Xinghong He                 Last update: 05 Aug 1988   *
!   Modified for ifort -i8 by A. Kramida (AK) Last update 22 Mar 2016  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE CORRE_C
      USE damp_C
      USE def_C
      USE grid_C
      USE int_C
      USE mpi_C
      USE orb_C
      USE orthct_C
      USE scf_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cofpotmpi_I
      USE defcor_I
      USE solve_I
      USE orthsc_I
!GG      USE matrixmpi_I
!GG      USE newcompi_I
      USE setlagmpi_I
      USE quad_I
      USE consis_I
      USE dampck_I
      USE dampor_I
      USE orthy_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J
      REAL(DOUBLE), INTENT(INOUT) :: DAMPMX
      LOGICAL  :: EOL, LSORT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: P2 = 2.0D-01
      REAL(DOUBLE), PARAMETER :: P005 = 5.0D-03
      REAL(DOUBLE), PARAMETER :: P0001 = 1.0D-04
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IPR, NPTS, INV, JP, NNP, I, NWWW
      INTEGER :: I_MPI, ndcof_max, ntot, i_last, iproc, ndcip, jproc
      INTEGER :: ifound, ind_buf, K
      REAL(DOUBLE) :: ED1, GAMAJ, ED2, WTAEV, DNORM, DNFAC, DEL1, DEL2, ODAMPJ
      LOGICAL :: FAIL, FIRST
      REAL(DOUBLE), DIMENSION(:), POINTER :: da_buffer
      INTEGER, DIMENSION(:), POINTER :: nda_buffer,ndcof_buffer
!-----------------------------------------------
!
!   C Froese Fischer's IPR and ED1 parameter
!
      DATA IPR/ 0/
      DATA ED1/ 0.D0/
      DATA FIRST/ .FALSE./
!
!
!-----------------------------------------------------------------------
!AK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
!      ISIZE = sizeof(NCF)
      I_MPI = MPI_INTEGER
!      if (ISIZE.EQ.8) I_MPI = MPI_INTEGER8
!
      GAMAJ = GAMA(J)
!
!   C Froese Fischer's parameters IPR, ED1, ED2 are set and
!   used in this routine and in DAMPCK
!
    1 CONTINUE
      ED2 = E(J)
!
!   Set up the exchange potential and arrays XU, XV as appropriate
!
!   Set coefficients for YPOT, XPOT, DACON
!   Compute direct potential, exchange potential
!   Add in Lagrange-multiplier contribution
!   Add in derivative-terms contribution
!
      NPTS = N
      CALL COFPOTmpi (EOL, J, NPTS)
!
!   Calculate deferred corrections
!
      CALL DEFCOR (J)
!
!   Solve the Dirac equation
!
!AK     NDCOF, NDA, and DA cannot be reduced by finding the max of NDCOF
!          and summing up DA from all nodes. The correct replacement is below
!cjb MPI_ALLREDUCE NDCOF
!      call MPI_ALLREDUCE(ndcof,ndcof_buffer,1,
!     :          I_MPI,MPI_MAX,MPI_COMM_WORLD,ierror)
!      if (ndcof_buffer .gt. ndcof) then
!        print *, 'improvmpi: ndcof, ndcof_buffer, myid=',
!     &                     ndcof, ndcof_buffer, myid
!        write (*,'(A35,4I10)') 'improv11: myid,NDCOF,ndcbfr=',
!     &    myid,NDCOF,ndcof_buffer
!        do indcof = ndcof+1, ndcof_buffer
!           da(indcof) = 0.0
!           nda(indcof) = 0
!        enddo
!        ndcof = ndcof_buffer
!      endif
      INV = 0
!cjb MPI_ALLREDUCE DA
!      if (ndcof.gt.0) then
!         call alloc(pda_buffer, ndcof, 8)
!         call alloc(pnda_buffer, ndcof, ISIZE)
!         call MPI_ALLREDUCE(da,da_buffer,ndcof,
!     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
!         da(1:ndcof) = da_buffer(1:ndcof)
!      nda(1:ndcof) = nda_buffer(1:ndcof)
!         call dalloc(pda_buffer)
!         call dalloc(pnda_buffer)
!      ENDIF
!AK      Consolidate NDA and DA from all nodes
      if (nprocs.gt.1) then
         call alloc(ndcof_buffer,nprocs,'ndcof_buffer','IMPROVmpi')
         call MPI_GATHER(ndcof, 1, I_MPI, ndcof_buffer,            &
              1, I_MPI, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_Bcast (ndcof_buffer, nprocs, I_MPI, 0,           &
                               MPI_COMM_WORLD, ierr)
         ndcof_max = 0
         do i = 1, nprocs
            if (ndcof_buffer(i).gt.ndcof_max) ndcof_max=ndcof_buffer(i)
         ENDDO
         if (ndcof_max.gt.0) then
            ntot = ndcof_max*nprocs
            call alloc(nda_buffer,ntot,'nda_buffer','IMPROVmpi')
            call MPI_GATHER(nda, ndcof_max, I_MPI, nda_buffer,     &
                 ndcof_max, I_MPI, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_Bcast (nda_buffer, ntot, I_MPI, 0,            &
                               MPI_COMM_WORLD, ierr)
            call alloc(da_buffer, ntot, 'da_buffer','IMPROVmpi')
            call MPI_GATHER(da, ndcof_max, MPI_DOUBLE_PRECISION,   &
              da_buffer, ndcof_max, MPI_DOUBLE_PRECISION, 0,       &
              MPI_COMM_WORLD, ierr)
            CALL MPI_Bcast (da_buffer, ntot, MPI_DOUBLE_PRECISION, &
                               0, MPI_COMM_WORLD, ierr)
            i_last = 0
            do iproc = 1, nprocs
               ndcip = ndcof_buffer(iproc)
               do jproc = 1, ndcip
                  ifound = 0
                  ind_buf = (iproc - 1)*ndcof_max + jproc
                  do k = 1, i_last
                     if (nda_buffer(ind_buf) .eq. NDA(k)) THEN
                        DA(k) = DA(k) + da_buffer(ind_buf)
                        ifound = k
                        exit
                     endif
                  enddo
                  if (ifound.eq.0) then
                     i_last = i_last + 1
                     do while (i_last .GT. NDDIM)
                        IF (NDDIM .GT. 0) THEN
                           NDDIM = 2*NDDIM
!GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 2)
                           CALL RALLOC(NDA,NDDIM,'NDA','IMPROVmpi')
                        ELSE
                          NDDIM = 64
!GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 1)
                           CALL ALLOC(NDA,NDDIM,'NDA','IMPROVmpi')
                        ENDIF
                     ENDdo
                     NDA(i_last) = nda_buffer(ind_buf)
                     DA(i_last) = da_buffer(ind_buf)
                  endif
               enddo
            enddo
            NDCOF = ndcof_max
            call dalloc(nda_buffer,'nda_buffer','IMPROVmpi')
            call dalloc(da_buffer,'da_buffer','IMPROVMPI')
         endif
         call dalloc(ndcof_buffer,'ndcof_buffer','IMPROVmpi')
      endif

      CALL SOLVE (J, FAIL, INV, JP, NNP)
!
!   Upon failure issue message; take corrective action if possible
!
      IF (FAIL) THEN
         IF (MYID == 0) WRITE (*, 300) NP(J), NH(J), METHOD(J)
         IF (METHOD(J) /= 2) THEN
            METHOD(J) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            CALL ORTHSC
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN
!               CALL MATRIXMPI
!               CALL NEWCOMPI (WTAEV)
!            ENDIF
            CALL SETLAGmpi (EOL)
            GO TO 1
         ELSE
            IF (MYID == 0) WRITE (*, 301)
            !CALL TIMER (0)
            STOP
         ENDIF
      ENDIF
!
!   Compute norm of radial function
!
      TA(1) = 0.D0
      TA(2:MTP0) = (P(2:MTP0)**2+Q(2:MTP0)**2)*RP(2:MTP0)
      MTP = MTP0

      CALL QUAD (DNORM)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      CALL CONSIS (J)
!
!   Normalize
!
      DNFAC = 1.D0/DSQRT(DNORM)
      P0 = P0*DNFAC
      P(:MTP0) = P(:MTP0)*DNFAC
      Q(:MTP0) = Q(:MTP0)*DNFAC
!
!   Check if different method should be used or if improvement
!   count should be reduced
!
      DEL1 = DABS(1.D0 - ED2/E(J))
      IF (METHOD(J) == 1) THEN
         DEL2 = DMAX1(DABS(1.D0 - DSQRT(DNORM)),DABS(DNFAC - 1.D0))
         IF (DEL1<P005 .AND. DEL2>P2) THEN
            METHOD(J) = 2
            GO TO 1
         ENDIF
      ELSE
         IF (DEL1<P0001 .AND. NSIC>1) NSIC = NSIC - 1
      ENDIF
!
!   Damp the orbital --- if not converged
!
      IF (SCNSTY(J) > ACCY) THEN
         CALL DAMPCK (IPR, J, ED1, ED2)
         ODAMPJ = DABS(ODAMP(J))
      ELSE
         ODAMPJ = 0.D0                           ! take the whole new orbital
      ENDIF
      CALL DAMPOR (J, INV, ODAMPJ)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      IF (ORTHST) THEN
         !CALL orthor (J, inv)
         NWWW = NW
         CALL ORTHY (NWWW, J, LSORT)
      ENDIF
!
!   Print details of iteration
!
      IF (MYID == 0)                                                &
         WRITE (*, 302) NP(J),NH(J),E(J),METHOD(J),PZ(J),SCNSTY(J), &
!cjb DNORM-1 -> SQRT(DNORM)-1
!cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
                        SQRT(DNORM)-1, ODAMPJ, JP, MF(J), INV, NNP
      DAMPMX = DMAX1(DAMPMX,DABS(ODAMPJ))

  300 FORMAT(/,' Failure; equation for orbital ',1I2,1A2,&
         ' could not be solved using method ',1I1)
  301 FORMAT(/,/,' ****** Error in SUBROUTINE IMPROV ******'/,&
         ' Convergence not obtained'/)
  302 FORMAT (1X,1I2,1A2,1P,1D16.7,1x,1I2,D11.3,1D10.2,1D10.2,&
              0P,F6.3,1x,1I5,1x,1I5,1x,1I2,1x,1I2)
      RETURN
      END SUBROUTINE IMPROVmpi
