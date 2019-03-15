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
      REAL(DOUBLE)          :: DENSTY
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
         WRITE (24, 301) CUTOFFTMP
         WRITE (24, 302) NCOEITMP
         WRITE (24, 303) NCOECTMP
         WRITE (24, 304) NCTEITMP
         WRITE (24, 305) NCTECTMP
         IF (LTRANS) THEN
            WRITE (24, 306) NTPITMP
            WRITE (24, 307) NMCBPTMP
            WRITE (24, 308) NCORETMP
         ENDIF
         IF (LVP) WRITE (24, 309) NVPITMP
         IF (LNMS) WRITE (24, 310) NKEITMP
         IF (LSMS) WRITE (24, 311) NVINTITMP
      ELSE
         WRITE (24, *) 'Restart mode --- no report on radial integrals'
      ENDIF                                      !(irestart .NE. 0) ! non-restart mode


!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------

      NELMNT_A = NELMNTTMP

      DENSTY = DBLE(NELMNTTMP)/DBLE((NCFTMP*(NCFTMP + 1))/2)

      EAV = EAV/DBLE(NCFTMP) + ELSTO

      WRITE (24, 312) NELMNTTMP
      WRITE (24, 313) DENSTY
      WRITE (24, *)
      WRITE (24, 300) EAV
      WRITE (24, *) EAV, ELSTO, NCFTMP

  300 FORMAT('Average energy = ',1P,D19.12,' Hartrees.')
  301 FORMAT('CUTOFF set to ',1P,D17.10)
  302 FORMAT('Dirac-Coulomb one-e radial integrals:',1I8)
  303 FORMAT('One-e angular integrals that exceed CUTOFF: ',1I8)
  304 FORMAT('Coulomb two-e radial integrals: ',1I8)
  305 FORMAT('Two-e angular integrals that exceed CUTOFF: ',1I8)
  306 FORMAT('Transverse two-e radial integrals: '/,6I8)
  307 FORMAT('MCBP coefficients that exceed CUTOFF: ',1I8)
  308 FORMAT('Core coefficients that exceed CUTOFF: ',1I8)
  309 FORMAT('Vacuum polarisation integrals: ',1I8)
  310 FORMAT('Kinetic energy integrals: ',1I8)
  311 FORMAT('Vinti integrals: ',1I8)
  312 FORMAT('Elements that exceed CUTOFF in the lower',&
         ' triangle of the H matrix: ',1I8)
  313 FORMAT('Density of the H(amiltonian) matrix: ',1P,D22.15)

      RETURN
      END SUBROUTINE GENMAT2
