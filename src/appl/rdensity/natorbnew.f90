!***********************************************************************
!                                                                      *
      SUBROUTINE NATORBNEW(NAME,DENS1VEC,DINT1VEC)
!                                                                      *
!   IF angular coefficients must be calculated                         *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the electron density and the natural orbitals       * 
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
!                        ITJPO, ONESCALAR                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Feb. 2012                                   *
!   Modified by J. Ekman   Nov. 2013                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW, NNNP
      USE orb_C
      USE debug_C,         ONLY: CUTOFF
      USE wave_C
      USe grid_C
      USE prnt_C
      USE eigv_C
      USE syma_C
      USE IOUNIT_C,     ONLY: ISTDI, ISTDE !SACHA
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE convrt_I
      USE onescalar_I
      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!      REAL(DOUBLE) :: dsyev
      EXTERNAL     :: dsyev
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (LEN = 24), INTENT(IN) :: NAME
      REAL(DOUBLE), DIMENSION(NVEC,NNNP), INTENT(OUT)     :: DENS1VEC
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NNNP), INTENT(IN) :: DINT1VEC
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*11  :: CNUM
      CHARACTER*20  :: ORBKAPPA(-10:10)                          
      CHARACTER*256 :: FILNAM                                
      REAL(DOUBLE), DIMENSION(:,:,:,:), ALLOCATABLE :: DENSMAT          
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: DENSMAT2            
      REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: DENSMAT3           
      REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: eig,work            
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NORB                  
      INTEGER, DIMENSION(:), ALLOCATABLE :: KFLAG                  
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL_S, IA_S
      REAL(DOUBLE), DIMENSION(100) :: DUMMY                     
      REAL(DOUBLE) :: SWSUM, CONTRI1, ELEMNT1,CONTRIB
      INTEGER, DIMENSION(100) :: IFIRST,NC                     
      INTEGER, DIMENSION(NNNW) :: ROTCHECK                    
      INTEGER :: ERROR, MAXK, MAXN, MFMAX, MINK, MINN, LPAR, LOC
      INTEGER :: K, KA, NCONTR, IR, IC, IK, IV, LAB, IOPAR, IA, IB,IAA,ISTORE
      INTEGER :: INCOR, I, J, L, NDIM, ITJPOC, LCNUM, inf, REORDERING,PURE_EGVC
!-----------------------------------------------
!
	WRITE(ISTDE,*)' How do you want to order your egvc ?'
	WRITE(ISTDE,*) '                1) By looking at the dominant component'
	WRITE(ISTDE,*) '                2) Following the decreasing order of the egvl'
        REORDERING=0
	DO WHILE (REORDERING .NE. 1 .AND. REORDERING .NE. 2) 
	 read(ISTDI,*) REORDERING
	ENDDO


      MAXK = maxval(NAK(1:NW))
      MINK = minval(NAK(1:NW))
      MAXN = maxval(NP(1:NW))
      MINN = minval(NP(1:NW))

      ORBKAPPA(-7) = '7i 8i ..'
      ORBKAPPA(-6) = '6h 7h ..'
      ORBKAPPA(-5) = '5g 6g ..'
      ORBKAPPA(-4) = '4f 5f ..'
      ORBKAPPA(-3) = '3d 4d ..'
      ORBKAPPA(-2) = '2p 3p ..'
      ORBKAPPA(-1) = '1s 2s ..'
      ORBKAPPA(1) = '2p- 3p- ...'
      ORBKAPPA(2) = '3d- 4d- ...'
      ORBKAPPA(3) = '4f- 5f- ...'
      ORBKAPPA(4) = '5g- 6g- ...'
      ORBKAPPA(5) = '6h- 7h- ...'
      ORBKAPPA(6) = '7i- 7i- ...'

      write(*,*) 'Maximum Kappa: ', MAXK
      write(*,*) 'Minimum Kappa: ', MINK
      write(*,*) 'Maximum N: ', MAXN
      write(*,*) 'Minimum N: ', MINN
!      write(*,*)
!      write(*,*) 'DIAG/NDIAG CSF1        CSF2        2J+1       Kappa        
!     :      nl             nl   Contribution to density matrix'
!      write(*,*) '------------------------------------------------------
!     :----------------------------------------------------'

      ALLOCATE( DENSMAT(NVEC,MINK:MAXK,MINN:MAXN,MINN:MAXN) )
      ALLOCATE( DENSMAT2(20,20) )
      ALLOCATE( DENSMAT3(20,20) )
      ALLOCATE( eig(MAXN) )
      ALLOCATE( NORB(MINK:MAXK,MINN:MAXN) )
      ALLOCATE( KFLAG(MINK:MAXK) )

      ROTCHECK(:) = 0             
      DENS1VEC(:,:) = 0.0D00     
      DENSMAT(:,:,:,:) = 0.0D00 

      KFLAG(:) = 0
      DO I=1,NW
         KFLAG(NAK(I)) = KFLAG(NAK(I)) + 1
      END DO

!     
!     Set the rank (zero) and parity (even) for the one-particle
!     coefficients
!     
      KA = 0
      IOPAR = 1
      INCOR = 1
!     
!     Allocate storage for the arrays in BUFFER
!     
      CALL ALCBUF (1)
!     
!     Sweep through the Hamiltonian matrix to determine the
!     sms parameter
!     
      DO 13 IC = 1,NCF
!     
!     Output IC on the screen to show how far the calculation has preceede
!     
         CALL CONVRT (IC,CNUM,LCNUM)
         if (mod(IC,100).eq.0) then
            PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
         end if
!     
         ITJPOC = ITJPO(IC)
         DO 12 IR = IC,NCF
!            write(*,*) 'IC: ', IC, '  IR: ', IR
!            write(*,*) '--------------------------------------------'
!     
!     Matrix elements are diagonal in J
!     
            IF (ITJPO(IR) .EQ. ITJPOC) THEN
!     
!     Initialise the accumulator
!     
               ELEMNT1 = 0.0D00
!     
!     Call the MCT package to compute T coefficients
!     
               CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
               IF (IA .NE. 0) THEN
                  IF (IA .EQ. IB) THEN
                     NCONTR = 0
                     DO 8 IA = 1,NW
                        IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                           NCONTR = NCONTR + 1
                           TSHELL_S(NCONTR) = TSHELL(IA)
                           IA_S(NCONTR) = IA
                           ELEMNT1 = TSHELL(IA)
                           DO 9 J = 1,NVEC
                              LOC = (J-1)*NCF
                             CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
                              ! FILL DENSITY MATRIX WITH DIAGONAL ELEMENTS
                              DENSMAT(J,NAK(IA),NP(IA),NP(IA)) =       &
                              DENSMAT(J,NAK(IA),NP(IA),NP(IA)) + CONTRI1
			      IF (IC .NE. IR) CONTRI1=CONTRI1*2.0D00 
			      DO 23 L = 2,NNNP                                    
                      		DENS1VEC(J,L) = DENS1VEC(J,L)+DINT1VEC(IA,IA,L)*& 
                                      CONTRI1                         
   23               	      CONTINUE 
 9                         CONTINUE
                        ENDIF
 8                   CONTINUE
                  ELSE          !IF IA NOT EQUAL TO IB
                     IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                        IF (NAK(IA).EQ.NAK(IB)) THEN
                           ELEMNT1 = TSHELL(1)
                           DO 10 J = 1,NVEC
                              LOC = (J-1)*NCF
                             CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
                              ! FILL DENSITY MATRIX WITH NON-DIAGONAL ELEMENTS
                              DENSMAT(J,NAK(IA),NP(IA),NP(IB)) =       &
                              DENSMAT(J,NAK(IA),NP(IA),NP(IB)) + CONTRI1
                              DENSMAT(J,NAK(IA),NP(IB),NP(IA)) =       &
                              DENSMAT(J,NAK(IA),NP(IA),NP(IB))
			   IF (IR .NE. IC) CONTRI1=CONTRI1*2.0D00
			   DO 21 L = 2,NNNP                                    
                      		DENS1VEC(J,L) = DENS1VEC(J,L)+DINT1VEC(IA,IB,L)*& 
                                      CONTRI1                        
   21               		CONTINUE
 10                        CONTINUE
                        ENDIF   ! END KAPPA(IA) EQUAL TO KAPPA(IB)
                     ENDIF      ! END CUTOFF
                  ENDIF         ! END IA NOT EQUAL TO IB
               ENDIF            ! END IA NOT EQUAL TO 0
            ENDIF  ! ENDS  IF (ITJPO(IR) .EQ. ITJPOC) THEN
 12      CONTINUE
 13   CONTINUE


      ! CONSTRUCT ARRAY THAT GIVES THE ORBIT ORDER AS A 
      ! FUNCTION OF KAPPA AND N. TO BE USED BELOW!
      DO I = 1,NW
         NORB(NAK(I),NP(I)) = I
      END DO
         
      ! LOOP OVER KAPPA FROM KAPPA = MINK TO KAPPA = MAXK
      DO IK = MINK,MAXK
         ! KAPPA = 0 DOES NOT EXIST
         IF (KFLAG(IK).GT.0) THEN
            ! DETERMINE MINN AND MAXN FOR KAPPA = IK
            MINN = 100
            MAXN = 0
            DO I = 1,NW
               IF(NAK(I).EQ.IK) THEN
                  IF(NP(I).LT.MINN) MINN = NP(I)
                  IF(NP(I).GT.MAXN) MAXN = NP(I)
               END IF
            END DO

            !IF (IK.LT.0) MINN = ABS(IK)
            !IF (IK.GT.0) MINN = IK + 1
            
            NDIM = MAXN - MINN + 1
            LPAR=NDIM*(3+NDIM/2)
            ALLOCATE(work(LPAR)) 
            write(*,*) 
            write(*,*) '-------------------------------------'
            write(*,'(a,i4,a,a)') 'KAPPA: ', IK, '  ORB: ', ORBKAPPA(IK)
            write(*,*) '-------------------------------------'
            write(*,*) 'NDIM: ', NDIM

            ! LOOP OVER STATES AND FILL DENSITY MATRIX DENSMAT2 FOR KAPPA = IK
            ! MATRICES ARE WEIGHTED WITH THEIR STATISTICAL WEIGHT: 2J+1 (IATJPO)
            SWSUM = 0.0d0
            DO IV = 1,NVEC
               DO IA = 1,NDIM
                  DO IB = 1,NDIM
                     IF(IV.EQ.1) THEN
                        DENSMAT2(IA,IB) = IATJPO(IV)*                  &
                             DENSMAT(IV,IK,IA+MINN-1,IB+MINN-1)
                     ELSE
                        DENSMAT2(IA,IB) = DENSMAT2(IA,IB) + IATJPO(IV)*&
                             DENSMAT(IV,IK,IA+MINN-1,IB+MINN-1)
                     END IF
                  END DO
               END DO
               SWSUM = SWSUM + IATJPO(IV)
            END DO

            !write(*,*) SWSUM

            ! OBTAIN THE WEIGHTED AVERAGE DENSITY MATRIX
            DENSMAT2 = DENSMAT2/SWSUM

            ! WRITE OUT DENSITY MATRIX FOR KAPPA = IK FOR 
            ! N = MINN TO N = MAXN
            write(*,*) 'density matrix:'
            DO IA = 1,NDIM
               write(*,'(i4,20f12.7)')IA+MINN-1,DENSMAT2(1:NDIM,IA)
            END DO
            write(*,*)
            
            ! CALL DSYEV FOR EIGENVALUES AND EIGENVECTORS
            ! NOTE THAT EIGENVALUES ARE GIVEN IN ASCENDING 
            ! ORDER AND EIGENVECTORS ARE GIVEN IN AN ORDER
            ! CORRESPONDING TO THE EIGENVALUES.
            ! EIGENVALUES ARE STORED IN ARRAY eig
            ! EIGENVECTORS ARE STORED IN 2D INPUT ARRAY DENSMAT2

            DENSMAT3(1:NDIM,1:NDIM) = DENSMAT2(1:NDIM,1:NDIM)
            call dsyev('V','U',NDIM,DENSMAT3(1:NDIM,1:NDIM),           &
                 NDIM,eig,work,LPAR,inf)
            
            ! WRITE OUT EIGENVALUES
            write(*,*) 'eigenvalues:'
            DO IA = 1,NDIM
               write(*,'(i4,20f12.7)') IA, eig(IA)
            END DO
            write(*,*)
            ! WRITE OUT EIGENVECTORS
            write(*,*) 'eigenvectors:'
	    DO IA=1,NDIM
		write(*,'(i4,20f12.7)') IA,DENSMAT3(1:NDIM,IA)
	    ENDDO
	    NC=0
	    IF (REORDERING .EQ. 1) THEN
	    write(*,*) 'eigenvectors after reordering according to dominance:'
            DO IA = 1,NDIM
		! WHEN IDENTICAL EIGENVALUES, E.G. MULTIPLE CLOSE SHELLS, WE NEED TO SORT THE EGVC
		! Should be use carefully !! Not working when two egvc have same dominance
		CONTRIB=0.d0
	       DO IB = 1, NDIM
		IF (ABS(DENSMAT3(IB,IA)) .EQ. 1.0) THEN
			NC(IB)=IA
			GO TO 50
                ELSEIF (ABS(DENSMAT3(IB,IA)) .GT. CONTRIB) THEN
                         CONTRIB=ABS(DENSMAT3(IB,IA))
                         ISTORE=IB
		ENDIF
	       ENDDO
		IF (NC(ISTORE) .NE. 0 ) THEN
			WRITE(*,*) ' '
			WRITE(*,*) 'Two EGVC with the same dominant component '
			WRITE(*,*) 'Rerun with the ordering of EGVC according to decreasing EGVL'
			WRITE(*,*) ' '
			STOP
		ENDIF	
		NC(ISTORE)=IA
 	50      continue
            END DO
	    DO IA = 1,NDIM
		write(*,'(2i4,20f12.7)') IA,NC(IA),DENSMAT3(1:NDIM,NC(IA))
	    ENDDO
            ELSEIF (REORDERING .EQ. 2) THEN
	    write(*,*) 'eigenvectors after reordering according to decreasing EGVL:'
            DO IA = 1,NDIM
		 !Should just be careful when some orbitals are closed (=pure in their egvc composition) and others are mixed
	         PURE_EGVC=0
		 DO IB = 1, NDIM
        	        IF (ABS(DENSMAT3(IB,IA)) .EQ. 1.0) THEN
                 	       NC(IB)=IA
                               PURE_EGVC=1
			ENDIF
		ENDDO
		IF (PURE_EGVC .EQ. 0) THEN
		NC(NDIM-IA+1)=IA
		ENDIF
	    ENDDO
	    DO IA = 1,NDIM
                write(*,'(i4,20f12.7)') IA,DENSMAT3(1:NDIM,NC(IA))
            ENDDO
	    ENDIF
            ! IF MATRIX DIMENSION > 1
            IF(NDIM.GT.1) THEN
               ! MATCH EIGENVECTORS WITH CORRESPONDING ORBIT (WITH PRINCIPAL QUANTUM NUMBER N)
               ! NC(IB=N-MINN+1) = IA (ROW IN EIGENVALUE MATRIX)
	       ! NEED SOME REORDERING BECAUSE THE DIAGONALIZATION SORTS THE EGVL IN INCREASING ORDER.
	       ! WE WANT TO SORT THEM ACCORDING TO 1) THE DOMINANT COMPONENT OF THE EGVC (=NOT FOLLOWING THE EGVL ORDERING)
	       !				   2) THE DECREASE ORDER OF EGVL
	       ! ONE SHOULD JUST BE CAREFUL TO CLOSED SHELLS WHERE THE EGVC IS EXACTLY 1.0. WHEN MULTIPLE CLOSED SHELLS ARE THERE, JUST SORT THEM ACCORDING TO DOMINANCE, 
	       ! SO THAT THE 1s REMAINS THE 1s IN THE NATURAL ORBITAL BASIS
               NC = 0
	       IF (REORDERING .EQ. 1) THEN
               DO IA = 1,NDIM
		  CONTRIB=0.d0
		  DO IB = 1,NDIM
		 	IF (ABS(DENSMAT3(IB,IA)) .EQ. 1.0) THEN
				NC(IB)=IA
				GO TO 51
                        ELSEIF (ABS(DENSMAT3(IB,IA)) .GT. CONTRIB) THEN
				CONTRIB=ABS(DENSMAT3(IB,IA))
				ISTORE=IB
			ENDIF
		  ENDDO
		  NC(ISTORE)=IA
       51         continue
               END DO
	       ELSEIF (REORDERING .EQ. 2) THEN
		DO IA = 1,NDIM
		  PURE_EGVC=0
                  DO IB = 1,NDIM
                        IF (ABS(DENSMAT3(IB,IA)) .EQ. 1.0) THEN
                                NC(IB)=IA
  				PURE_EGVC=1
                        ENDIF
                  ENDDO
		IF (PURE_EGVC .EQ. 0) THEN
                	NC(NDIM-IA+1)=IA
                ENDIF
               END DO
	       ENDIF
               ! IDENTIFY ORBITS - IFIRST(IA) = NORB(IK,MINN+IA-1) AND
               ! DETERMINE MAXIMUM NUMBER OF RADIAL POINTS OCCURING
               ! IN ORBITS CORRESPONDING TO THE EIGENVECTORS
               MFMAX = 0
               DO IA = 1,NDIM
                  IFIRST(IA) = NORB(IK,MINN+IA-1)
                  ROTCHECK(IFIRST(IA)) = 1
                  IF(MF(IFIRST(IA)).GT.MFMAX) THEN
                     MFMAX = MF(IFIRST(IA))
                  END IF
               END DO
               ! LOOP OVER RADIAL POINTS
               DO K = 1,MFMAX
                  ! STORE LARGE COMPONENT (P) WAVEFUNCTION VALUE OF POINT K
                  ! IN ARRAY DUMMY FOR THE NDIM ORBITALS
                  DO IA = 1,NDIM
                     DUMMY(IA) = PF(K,IFIRST(IA))
                  END DO
                  ! DETERMINE THE NATURAL P WAVEFUNCTION VALUE OF 
                  ! POINT K FOR THE NDIM ORBITALS USING THE 
                  ! EIGENVECTORS
                  DO IA = 1,NDIM
                     DO IB = 1,NDIM
                        IF(IB.EQ.1) THEN
                           PF(K,IFIRST(IA)) =                      &
                                DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                        ELSE
                           PF(K,IFIRST(IA))=PF(K,IFIRST(IA)) + &
                                DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                        ENDIF
                     END DO
                  END DO
                  ! STORE SMALL COMPONENT (Q) WAVEFUNCTION VALUE OF POINT K
                  ! IN ARRAY DUMMY FOR THE NDIM ORBITALS
                  DO IA = 1,NDIM
                     DUMMY(IA) = QF(K,IFIRST(IA))
                  END DO
                  ! DETERMINE THE NATURAL Q WAVEFUNCTION VALUE OF 
                  ! POINT K FOR THE NDIM ORBITALS USING THE 
                  ! EIGENVECTORS
                  DO IA = 1,NDIM
                     DO IB = 1,NDIM
                        IF(IB.EQ.1) THEN
                           QF(K,IFIRST(IA)) =                      &
                                DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                        ELSE
                           QF(K,IFIRST(IA))=QF(K,IFIRST(IA)) + &
                                DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                        ENDIF
                     END DO
                  END DO
               END DO                  ! LOOP OVER RADIAL POINTS ENDS
               ! STORE PZ VALUE IN ARRAY DUMMY FOR THE NDIM ORBITALS
               DO IA = 1,NDIM
                  DUMMY(IA) = PZ(IFIRST(IA))
               ! DETERMINE PZ VALUES OF THE NDIM NATURAL ORBITS
               END DO
               DO IA = 1,NDIM
                  DO IB = 1,NDIM
                     IF(IB.EQ.1) THEN
                        PZ(IFIRST(IA)) =                           &
                             DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                     ELSE
                        PZ(IFIRST(IA)) = PZ(IFIRST(IA)) +      &
                             DENSMAT3(IB,NC(IA))*DUMMY(IB) 
                     ENDIF
                  END DO
                  ! IF DERIVATIVE AT ORIGIN IS NEGATIVE CHANGE SIGN OF
                  ! PF, QF AND PZ
                  IF(PZ(IFIRST(NC(IA))).LT.0.0) THEN
                     PF(:,IFIRST(NC(IA))) =                            &
                          -PF(:,IFIRST(NC(IA)))
                     QF(:,IFIRST(NC(IA))) =                            &
                          -QF(:,IFIRST(NC(IA))) 
                     PZ(IFIRST(NC(IA))) = -PZ(IFIRST(NC(IA)))
                  END IF
               END DO
               ! SET THE MAXIMUM NUMBER OF RADIAL POINTS TO MFMAX
               ! FOR ALL ORBITALS   
               DO IA = 1,NDIM
                  MF(IFIRST(IA)) = MFMAX
               END DO
            ENDIF                      ! IF(NDIM.GT.1) ENDS
            DEALLOCATE(WORK)
         ENDIF                         ! IF (ABS(IK).GT.0) ENDS
      END DO                           ! LOOP OVER KAPPA ENDS

      ! WRITE NATURAL ORBITALS TO FILE <NAME>.nw
      write(*,*)
      write(*,*) '--------------------------------------'
      write(*,*) 'Natural orbits are written to file'
      write(*,*) ' nl     KAPPA     PZ           ROTATED'
      write(*,*) '--------------------------------------'
      K=INDEX(NAME,' ')
      FILNAM = NAME(1:K-1)//'.nw'
      OPEN(36,FILE=FILNAM,FORM='UNFORMATTED',STATUS = 'UNKNOWN')
      WRITE(36) 'G92RWF'
      DO I = 1,NW
         WRITE(*,'(i3,a,i8,f15.7,i8)') NP(I),NH(I),NAK(I),PZ(I),       &
              ROTCHECK(I)
         WRITE(36) NP(I),NAK(I),E(I),MF(I)
         WRITE(36) PZ(I),(PF(J,I),J = 1,MF(I)),(QF(J,I),J = 1,MF(I))
         WRITE(36) (R(J),J = 1,MF(I))
      ENDDO
      CLOSE(36)

!     
!     Deallocate storage for the arrays in BUFFER
      CALL ALCBUF (3)

      DEALLOCATE(DENSMAT, STAT = ERROR)
      DEALLOCATE(DENSMAT2, STAT = ERROR)
      DEALLOCATE(DENSMAT3, STAT = ERROR)

      RETURN
      END SUBROUTINE NATORBNEW
