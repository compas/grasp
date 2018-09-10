!***********************************************************************
!                                                                      *
      SUBROUTINE GENMAT2(IRESTART, NELMNT_A, ELSTO) 
!
!   Get eav and do writings to the summary file .csum
!   The mpi version (genmat2mpi) also gets nelmnt_a and elsto
!
! Xinghong He 98-06-15
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE, LONG
      USE hmat_C  !! rather than setham_to_genmat2
      USE decide_C
      USE eigv_C
      USE mpi_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)        :: IRESTART 
      INTEGER(LONG), INTENT(OUT) :: NELMNT_A 
      REAL(DOUBLE), INTENT(IN)   :: ELSTO 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(6) :: NTPI_A
      INTEGER               :: NCOEI_A, NCTEI_A, NCTEC_A, NMCBP_A,    &
                               NCORE_A, NVPI_A, NKEI_A, NVINTI_A,     &
                               NCOEC_A
      REAL(DOUBLE)          :: DENSTY, EAV_A
!-----------------------------------------------
!-----------------------------------------------------------------------
! Transfer parameters around and print out on node-0.
!  .Temporary variables (ending with "_a") are used to store the
!   sums over nodes;
!  .Outputs to stream 24 are done on node-0;
! The following parameters are accumulated in setham which will not
! contain the correct values in restart mode. And thus is skipped.
!-----------------------------------------------------------------------
      IF (IRESTART /= 0) THEN                    ! non-restart mode 
!         WRITE (24, 301) CUTOFFTMP 
!         WRITE (24, 302) NCOEITMP 
!         WRITE (24, 303) NCOECTMP 
!         WRITE (24, 304) NCTEITMP 
!         WRITE (24, 305) NCTECTMP 
!         IF (LTRANS) THEN 
!            WRITE (24, 306) NTPITMP 
!            WRITE (24, 307) NMCBPTMP 
!            WRITE (24, 308) NCORETMP 
!         ENDIF 
!         IF (LVP) WRITE (24, 309) NVPITMP 
!         IF (LNMS) WRITE (24, 310) NKEITMP 
!         IF (LSMS) WRITE (24, 311) NVINTITMP 
!      ELSE 
!         WRITE (24, *) 'Restart mode --- no report on radial integrals' 
      CALL MPI_Reduce (NCOEItmp, NCOEI_a, 1, MPI_INTEGER, MPI_SUM, 0,  &
                                       MPI_COMM_WORLD, ierr)
      CALL MPI_Reduce (NCOECtmp, NCOEC_a, 1, MPI_INTEGER, MPI_SUM, 0,  &
                                       MPI_COMM_WORLD, ierr)
      CALL MPI_Reduce (NCTEItmp, NCTEI_a, 1, MPI_INTEGER, MPI_SUM, 0,  &
                                       MPI_COMM_WORLD, ierr)
      CALL MPI_Reduce (NCTECtmp, NCTEC_a, 1, MPI_INTEGER, MPI_SUM, 0,  &
                                       MPI_COMM_WORLD, ierr)

      IF (LTRANS) THEN
         CALL MPI_Reduce (NTPItmp, NTPI_a, 6, MPI_INTEGER, MPI_SUM, 0, &
                                       MPI_COMM_WORLD, ierr)
         CALL MPI_Reduce (NMCBPtmp, NMCBP_a, 1, MPI_INTEGER, MPI_SUM,0,&
                                       MPI_COMM_WORLD, ierr)
         CALL MPI_Reduce (NCOREtmp, NCORE_a, 1, MPI_INTEGER,MPI_SUM,0, &
                                       MPI_COMM_WORLD, ierr)
      ENDIF

      IF (LVP)                                                         &
         CALL MPI_Reduce (NVPItmp, NVPI_a, 1, MPI_INTEGER, MPI_SUM, 0, &
                                       MPI_COMM_WORLD, ierr)

      IF (LNMS)                                                        &
         CALL MPI_Reduce (NKEItmp, NKEI_a, 1, MPI_INTEGER, MPI_SUM, 0, &
                                       MPI_COMM_WORLD, ierr)

      IF (LSMS)                                                        &
         CALL MPI_Reduce (NVINTItmp, NVINTI_a, 1, MPI_INTEGER, MPI_SUM,&
                                       0, MPI_COMM_WORLD, ierr)

      IF (myid .EQ. 0) THEN
         WRITE (24,301) CUTOFFtmp
         WRITE (24,302) NCOEI_a
         WRITE (24,303) NCOEC_a
         WRITE (24,304) NCTEI_a
         WRITE (24,305) NCTEC_a
         IF (LTRANS) THEN
            WRITE (24,306) NTPI_a
            WRITE (24,307) NMCBP_a
            WRITE (24,308) NCORE_a
         ENDIF
         IF (LVP) WRITE (24,309) NVPI_a
         IF (LNMS) WRITE (24,310) NKEI_a
         IF (LSMS) WRITE (24,311) NVINTI_a
      ENDIF
            ELSE
      IF (myid .EQ. 0) THEN
         WRITE (24,*) 'Restart mode --- no report on radial integrals'
      ENDIF
      ENDIF                                      !(irestart .NE. 0) ! non-restart mode 

 
 
!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------
 
      CALL MPI_Allreduce (NELMNTtmp, NELMNT_a, 1, MPI_INTEGER8,  &
!AK:      MPI_INTEGER8 required to match the sizes of NELMNT_a and NELMNTtmp
!CFF &    MPI_SUM, MPI_COMM_WORLD, ierr) ! want MAX not SUM
     &    MPI_MAX, MPI_COMM_WORLD, ierr)
!cjb  NELMNT_a decides between sparse and dense
!cjb  see comments in maneig.f90

      CALL MPI_Bcast (ELSTO, 1, MPI_DOUBLE_PRECISION, 0,         &
                            MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce (EAV, EAV_a, 1, MPI_DOUBLE_PRECISION,   &
                              MPI_SUM, MPI_COMM_WORLD, ierr)
      EAV = EAV_a / DBLE (ncftmp) + ELSTO

      IF (myid .EQ. 0) THEN
         WRITE (24,312) NELMNT_a
         WRITE (24, *)
         WRITE (24, 300) eav
      ENDIF 
 
  300 FORMAT('Average energy = ',1P,D19.12,' Hartrees.') 
  301 FORMAT('CUTOFF set to ',1P,D17.10) 
  302 FORMAT('Dirac-Coulomb one-e radial integrals:',1I8) 
  303 FORMAT('One-e angular integrals that exceed CUTOFF: ',1I8) 
  304 FORMAT('Coulomb two-e radial integrals: ',1I8) 
  305 FORMAT('Two-e angular integrals that exceed CUTOFF: ',1I11) 
  306 FORMAT('Transverse two-e radial integrals: '/,6I8) 
!cjb 1I8 -> 1I16
! 307 FORMAT('MCBP coefficients that exceed CUTOFF: ',1I8)
  307 FORMAT('MCBP coefficients that exceed CUTOFF: ',1I16) 
  308 FORMAT('Core coefficients that exceed CUTOFF: ',1I8) 
  309 FORMAT('Vacuum polarisation integrals: ',1I8) 
  310 FORMAT('Kinetic energy integrals: ',1I8) 
  311 FORMAT('Vinti integrals: ',1I8) 
  312 FORMAT('Elements that exceed CUTOFF in the lower',&
         ' triangle of the H matrix: ',1I11) 
  313 FORMAT('Density of the H(amiltonian) matrix: ',1P,D22.15) 
 
      RETURN  
      END SUBROUTINE GENMAT2 
