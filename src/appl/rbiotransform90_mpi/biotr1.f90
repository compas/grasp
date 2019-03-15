!***********************************************************************
!                                                                      *
      SUBROUTINE BIOTR1(PI, QI, NLI, NINSHLI, PF, QF, NLF, NINSHLF, &
             NGRID, MXL, SCR, LSCR, NTESTG, CISHL, CICI, CFSHL, CFCI)
!                                                                      *
! Generate Matrices for rotating radial functions and                  *
! for counter rotating CI coefficients                                 *
!                                                                      *
! This subroutine is a modification of the corresponding MCHF          *
! routine                                                              *
!                                                                      *
! Per Jonsson,     Lund June 1996                                      *
!                                                                      *
! This subroutine has been modified to support transformations         *
! in cases like 3s2p1d 3s2p                                            *
!                                                                      *
! Per Jonsson,     Lund Feb  1997                                      *
!                                                                      *
! =====                                                                *
! Input                                                                *
! =====                                                                *
!                                                                      *
! PI,QI   : Input radial shells for initial state                      *
!         : (large and small component)                                *
! NLI     : Number of shells per kappa for initial state               *
! NINSHLI : Number of inactive shells  per kappa for initial state     *
!                                                                      *
! PF,QF   : Input radial shells for final   state                      *
!         : (large and small component)                                *
! NLF     : Number of shells per kappa for final   state               *
! NINSHLF : Number of inactive shells  per kappa for final state       *
!                                                                      *
! NGRID   : Number of gridpoints                                       *
! MXL     : Number of kappa quantum numbers                            *
! SCR     :                                                            *
! LSCR    : Total length of scratch space.                             *
! NTESTG  : Global PRINT flag : = 0 => complete silence                *
!                               = 1 => test and PRINT overlap matrices,*
!                                      PRINT header and  t matrix      *
!                            .gt. 1 => Hope you now what you are doing *
!                                                                      *
! ======                                                               *
! Output                                                               *
! ======                                                               *
!                                                                      *
! CISHL : Rotation matrix for initial state shells                     *
! CICI  : Rotation matrix for initial state CI coefficients            *
! CFSHL : Rotation matrix for final state shells                       *
! CFCI  : Rotation matrix for final state CI coefficients              *
! PI,QI : Initial shells in biorthogonal basis                         *
!         biorthogonal basis                                           *
! PF,QF : Final   shells in biorthogonal basis                         *
!                                                                      *
!                                                                      *
! Inactive Shells : The functions of the inactive shells must be       *
!                   be supplied to the program, and NLI,NLF            *
!                   must refer to the total number of                  *
!                   occupied shells ( inactive+active)                 *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ifnmnx_I
      USE ielsum_I
!      USE gets_I
!      USE wrtmat_I
!      USE copvec_I
!      USE invmat_I
!      USE ulla_I
!      USE trpmat_I
!      USE matml4_I
!      USE scalve_I
!      USE setvec_I
!      USE pamtmt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NGRID,  MXL
      INTEGER, INTENT(IN) :: LSCR,  NTESTG
      INTEGER, DIMENSION(MXL)  :: NLI, NINSHLI, NLF, NINSHLF
!     The following have dimensions (NGRID:*)
!GG      REAL(DOUBLE), DIMENSION(:,:), pointer :: PI, QI, PF, QF
      REAL(DOUBLE), DIMENSION(NGRID,*) :: PI, QI, PF, QF
!     The following have fixed dimensions (NGRID, NLMAX)
      REAL(DOUBLE), DIMENSION(*)  :: CISHL, CFSHL
!     The work array
!     REAL(DOUBLE), DIMENSION(LWORK1) :: work
!     Th following have fixed dimensions (20*NLMAX*NLMAX)
      REAL(DOUBLE), DIMENSION(*) :: cici, cfci
      REAL(DOUBLE), DIMENSION(*) :: SCR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, ILI, ILF, NLIMX, NLFMX, NLIFMX, NLTI, &
                 NLTF, KFREE, KLPI, KLPF, KLPIF, KLSTOT, KLSIF, KLSIFI, &
                 KLCI, KLCF, KLSCR, L, IIOFF, IFOFF, NI, NF, NIFMN, III,&
                 JJJ, KLPMX, KLPMN, NMX, NMN, KLCMX, KLCMN, NDIFF, J, &
                 KLTI, I, KLTF
      REAL(DOUBLE) :: TII, TIII
!-----------------------------------------------
!
      NTESTL = 0
      NTEST = MAX(NTESTL,NTESTG)
      IF (NTEST >= 1) THEN
         WRITE (6, *)
         WRITE (6, *) '                   *************************'
         WRITE (6, *) '                   *   Entering BIOTR1     *'
         WRITE (6, *) '                   *************************'
         WRITE (6, *)
      ENDIF
!
!. Scratch should at least be of length
!
      ILI = 1
      ILF = 1
!
!. Largest number of shells of a given symmetry
!
      NLIMX = IFNMNX(NLI,MXL,1)
      NLFMX = IFNMNX(NLF,MXL,1)
      NLIFMX = MAX(NLIMX,NLFMX)
      IF (NTEST >= 10) WRITE (6, *) ' NLIMX,NLFMX NLIFMX ', NLIMX, NLFMX, &
         NLIFMX
!
!. Total numner of shells
!
      NLTI = IELSUM(NLI,MXL)
      NLTF = IELSUM(NLF,MXL)
      IF (NTEST >= 10) WRITE (6, *) ' NLTI NLTF', NLTI, NLTF
!
! Scratch space for orbital rotations
      KFREE = 1
!
      KLPI = KFREE
      KFREE = KFREE + NLIMX*NGRID
!      WRITE (6, *) ' In biotrn: KLPI    = ', KLPI
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLPF = KFREE
      KFREE = KFREE + NLFMX*NGRID
!      WRITE (6, *) ' In biotrn: KLPF    = ', KLPF
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLPIF = KFREE
      KFREE = KFREE + NLIFMX*NGRID
!      WRITE (6, *) ' In biotrn: KLPIF   = ', KLPIF
!      WRITE (6, *) '            KFREE   = ', KFREE
!
!. Total overlap matrix
!
      KLSTOT = KFREE
      KFREE = KFREE + NLTI*NLTF
!      WRITE (6, *) ' In biotrn: KLSTOT  = ', KLSTOT
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLSIF = KFREE
      KFREE = KFREE + NLIFMX**2
!      WRITE (6, *) ' In biotrn: KLSIF   = ', KLSIF
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLSIFI = KFREE
      KFREE = KFREE + NLIFMX**2
!      WRITE (6, *) ' In biotrn: KLSIFI  = ', KLSIFI
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLCI = KFREE
      KFREE = KFREE + NLIFMX**2
!      WRITE (6, *) ' In biotrn: KLCI    = ', KLCI
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLCF = KFREE
      KFREE = KFREE + NLIFMX**2
!      WRITE (6, *) ' In biotrn: KLCF    = ', KLCF
!      WRITE (6, *) '            KFREE   = ', KFREE
!
      KLSCR = KFREE
      KFREE = KFREE + NLIFMX**2 + NLIFMX*(NLIFMX + 1)
!      WRITE (6, *) ' In biotrn: KLSCR   = ', KLSCR
!      WRITE (6, *) '            KFREE   = ', KFREE
!      WRITE (6, *) '         =>  FREE =     ', KFREE
!
!. Check length of scratch
!
      IF (LSCR <= KFREE - 1) THEN
       WRITE (6,*)' BIOTR1 in trouble ! '
       WRITE (6,*)' Increase dimension of scratch before call to BIOTR1'
       WRITE (6,*)' Current and required length (LSCR,KFREE-1)', LSCR, &
            KFREE - 1
       STOP 'Increase LWORK before call to BIOTR1'
       ENDIF
!
!. Obtain overlap matrix
!
      CALL GETS (SCR(KLSTOT), NLTI,NLTF)

      DO L = 1, MXL
         IF (NTEST >= 5) THEN
            WRITE (6, *) '   L = ', L
            WRITE (6, *) '       Orbital rotation...'
         ENDIF
!
!. Offsets for given L in shell matrices
!
         IF (L == 1) THEN
            IIOFF = 1
            IFOFF = 1
         ELSE
            IIOFF = IIOFF + NLI(L-1)**2
            IFOFF = IFOFF + NLF(L-1)**2
         ENDIF
         IF (NTEST >= 1) THEN
           WRITE (6, *)
           WRITE (6, *) &
           ' BIOTRN : Information on transformations of shells with L =',L
           WRITE (6, *)
         ENDIF
!
! =========================================================
! 1 : Obtain Biorthogonal forms of initial and final shells
! =========================================================
!
!. The overlap matrix can be written
!
!     * * * * * * *
!     *       *   *
!     *   S   * X *
!     *       *   *
!     * * * * * * *
!
! where S is the quadratic subblock.
! The basis functions corresponding to the quadratic
! subblock are made biorthogonal with an UL decomposition.
! The remaining basis
! functions are made biorthogonal by choosing  the
! transformation
!            * * *
!            *   *
!            * Y *
!            *   *
!            * * *
!            * 1 *
!            * * *
!
!
! With Y = -S-1*X
!
         NI = NLI(L)
         NF = NLF(L)
!
!.1.1 : Obtain shells of given L in proper order for
!       biorthogonal treatment
!
         IF (L == 1) THEN
            ILI = 1
            ILF = 1
         ELSE
            ILI = ILI + NLI(L-1)
            ILF = ILF + NLF(L-1)
         ENDIF

!
! 1.2 obtain biorthogonal of the first min(ni,nf) shells
!
         NIFMN = MIN(NI,NF)
!
!. Overlap matrix SIF = Integral (PI(I)*PF(J))
!
!ww Per change to support cases like 3s2p1d 3s2p. The belonging endif
!ww Is just at the end

         IF (NIFMN <= 0) CYCLE

         DO III = 1, NIFMN
            DO JJJ = 1, NIFMN
               SCR(KLSIF+(JJJ-1)*NIFMN+III-1) = SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI&
                  +III+ILI-1)
            END DO
         END DO
         IF (NTEST >= 15) THEN
            WRITE (6, *) ' Overlap matrix '
            CALL WRTMAT (SCR(KLSIF), NIFMN, NIFMN, NIFMN, NIFMN)
         ENDIF
!
! Obtain upper triangular CI and CF so CI(T) S CF = 1
! or CF CI(T) = S-1, which corresponds to an UL decomposition
!
!. Invert S
!
         CALL COPVEC (SCR(KLSIF), SCR(KLSIFI), NIFMN**2)
         CALL INVMAT (SCR(KLSIFI), SCR(KLCI), NIFMN, NIFMN)

!. UL decompose
         CALL COPVEC (SCR(KLSIFI), SCR(KLSIF), NIFMN**2)
         CALL ULLA (SCR(KLSIF), SCR(KLCF), SCR(KLCI), NIFMN, SCR(KLSCR))
         CALL TRPMAT (SCR(KLCI), NIFMN, NIFMN, SCR(KLSCR))
         CALL COPVEC (SCR(KLSCR), SCR(KLCI), NIFMN**2)
!
!. The transformation matrix between the first NIFMX
!. shells is now known, biorthogonalize remaining orbitals
!
         IF (NI/=NF .AND. NI/=0 .AND. NF/=0) THEN
            IF (NI > NF) THEN
               KLPMX = KLPI
               KLPMN = KLPF
               NMX = NI
               NMN = NF
               KLCMX = KLCI
               KLCMN = KLCF
            ELSE
               KLPMX = KLPF
               KLPMN = KLPI
               NMX = NF
               NMN = NI
               KLCMX = KLCF
               KLCMN = KLCI
            ENDIF
            NDIFF = NMX - NMN
!
! Y = -S-1 * X
!. overlap X between remaining orbitals and the other set
!
            IF (NI > NF) THEN
!
! I columns F rows
!
               DO III = NMN + 1, NMX
                  DO JJJ = 1, NF
                     SCR(KLSIF+(III-NMN-1)*NF+JJJ-1) = SCR(KLSTOT-1+(JJJ+ILF-1-&
                        1)*NLTI+III+ILI-1)
                  END DO
               END DO
            ELSE IF (NF > NI) THEN
! F columns I rows
               DO JJJ = NMN + 1, NMX
                  DO III = 1, NI
                     SCR(KLSIF+(JJJ-NMN-1)*NI+III-1) = SCR(KLSTOT-1+(JJJ+ILF-1-&
                        1)*NLTI+III+ILI-1)
                  END DO
               END DO
            ENDIF
!
            IF (NI > NF) THEN
               CALL TRPMAT (SCR(KLSIFI), NMN, NMN, SCR(KLSCR))
               CALL COPVEC (SCR(KLSCR), SCR(KLSIFI), NMN**2)
            ENDIF
            CALL MATML4 (SCR(KLSCR), SCR(KLSIFI), SCR(KLSIF), NMN, NDIFF, NMN, &
               NMN, NMN, NDIFF, 0)
            CALL SCALVE (SCR(KLSCR), -1.0D0, NMN*NDIFF)
            CALL COPVEC (SCR(KLSCR), SCR(KLSIF), NMN*NDIFF)
!
! Construct complete CMX
!
            CALL SETVEC (SCR(KLSCR), 0.0D0, NMX**2)
            DO J = 1, NMX
               IF (J <= NIFMN) THEN
                  CALL COPVEC (SCR(KLCMX+(J-1)*NIFMN), SCR(KLSCR+(J-1)*NMX), &
                     NMN)
               ELSE
                  CALL COPVEC (SCR(KLSIF+(J-NMN-1)*NMN), SCR(KLSCR+(J-1)*NMX), &
                     NMN)
                  SCR(KLSCR-1+(J-1)*NMX+J) = 1.0D0
               ENDIF
            END DO
!
            CALL COPVEC (SCR(KLSCR), SCR(KLCMX), NMX**2)
         ENDIF
!ww Pertest
!         ENDIF
!
!. The two upper triangular matrices CI and CF are now known
!. Transfer to permanent arrays
!
         CALL COPVEC (SCR(KLCI), CISHL(IIOFF), NI**2)
         CALL COPVEC (SCR(KLCF), CFSHL(IFOFF), NF**2)
!
!. Rotate the large component of the shells
!
         CALL COPVEC (PI(1,ILI), SCR(KLPI), NI*NGRID)
         CALL COPVEC (PF(1,ILF), SCR(KLPF), NF*NGRID)
!
!         WRITE (*, *) 'Transformation matrices initial'
!         CALL WRTMAT (SCR(KLCI), NI, NI, NI, NI)
         CALL MATML4 (SCR(KLPIF), SCR(KLPI), SCR(KLCI), NGRID, NI, NGRID, NI, &
            NI, NI, 0)
         CALL COPVEC (SCR(KLPIF), PI(1,ILI), NI*NGRID)
!         WRITE (*, *) 'Transformation matrices final'
!         CALL WRTMAT (SCR(KLCF), NF, NF, NF, NF)
         CALL MATML4 (SCR(KLPIF), SCR(KLPF), SCR(KLCF), NGRID, NF, NGRID, NF, &
            NF, NF, 0)
         CALL COPVEC (SCR(KLPIF), PF(1,ILF), NF*NGRID)
!
!. Rotate the small component of the shells
!
         CALL COPVEC (QI(1,ILI), SCR(KLPI), NI*NGRID)
         CALL COPVEC (QF(1,ILF), SCR(KLPF), NF*NGRID)
!
         CALL MATML4 (SCR(KLPIF), SCR(KLPI), SCR(KLCI), NGRID, NI, NGRID, NI, &
            NI, NI, 0)
         CALL COPVEC (SCR(KLPIF), QI(1,ILI), NI*NGRID)
         CALL MATML4 (SCR(KLPIF), SCR(KLPF), SCR(KLCF), NGRID, NF, NGRID, NF, &
            NF, NF, 0)
         CALL COPVEC (SCR(KLPIF), QF(1,ILF), NF*NGRID)
!
         IF (NTEST >= 1) THEN
            WRITE (6, *) ' Test of overlap of biorthonormal functions'
! F columns I rows
            DO JJJ = 1, NF
               DO III = 1, NI
                  SCR(KLSIF+(JJJ-1)*NI+III-1) = SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI&
                     +III+ILI-1)
               END DO
            END DO
            CALL MATML4 (SCR(KLSCR), SCR(KLCI), SCR(KLSIF), NI, NF, NI, NI, NI&
               , NF, 1)
            CALL MATML4 (SCR(KLSIF), SCR(KLSCR), SCR(KLCF), NI, NF, NI, NF, NF&
               , NF, 0)
            WRITE (6, *) &
               ' new overlap matrix ( should be 1 on diag, 0 elsewhere )'
            CALL WRTMAT (SCR(KLSIF), NI, NF, NI, NF)
         ENDIF

         IF (NTEST >= 1) THEN
            WRITE (6, *)
            WRITE (6, *) ' Orbital Rotation matrix for I state'
            CALL WRTMAT (CISHL(IIOFF), NI, NI, NI, NI)
            WRITE (6, *) ' Orbital Rotation matrix for F state'
            CALL WRTMAT (CFSHL(IFOFF), NF, NF, NF, NF)
            WRITE (6, *)
         ENDIF
!
!. Matrix for counterrotation of CI coefficients, initial state
!
         KLTI = KLSIF
         CALL PAMTMT (SCR(KLCI), SCR(KLTI), SCR(KLSCR), NI)
         DO I = 1, NI
            TII = SCR(KLTI-1+(I-1)*NI+I)
            TIII = 1.0D0/TII
            CALL SCALVE (SCR(KLTI+(I-1)*NI), TIII, I - 1)
         END DO
         CALL COPVEC (SCR(KLTI), CICI(IIOFF), NI*NI)
!
!. Matrix for counterrotation of CI coefficients, Final state
!
         KLTF = KLSIF
         CALL PAMTMT (SCR(KLCF), SCR(KLTF), SCR(KLSCR), NF)
         DO I = 1, NF
            TII = SCR(KLTF-1+(I-1)*NF+I)
            TIII = 1.0D0/TII
            CALL SCALVE (SCR(KLTF+(I-1)*NF), TIII, I - 1)
         END DO
         CALL COPVEC (SCR(KLTF), CFCI(IFOFF), NF*NF)
         IF (NTEST < 1) CYCLE
         WRITE (6, *)
         WRITE (6, *) ' CI-Rotation matrix for I state'
         CALL WRTMAT (CICI(IIOFF), NI, NI, NI, NI)
         WRITE (6, *) ' CI-Rotation matrix for F state'
         CALL WRTMAT (CFCI(IFOFF), NF, NF, NF, NF)
         WRITE (6, *)

      END DO
!
      RETURN
      END SUBROUTINE BIOTR1
