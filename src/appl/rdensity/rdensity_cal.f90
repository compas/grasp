!***********************************************************************
!                                                                      *
      SUBROUTINE RDENSITY_CAL (NAME)
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the MS parameters, the electron density at the     *
!   origin and radial expectation values.                              *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, GETYN                          *
!               [SMS92]: RINTDENS, VINTI                               *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Feb. 2011                                   *
!   Modified by J. Ekman   Jan. 2014                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE memory_man
      USE debug_C
      USE decide_c
      USE def_C
      USE grid_C
      USE npar_C
      USE prnt_C
      USE syma_C
      USE orb_C
      USE teilst_C
      USE buffer_C
      USE rdensity_C
      USE jlabl_C, LABJ=> JLBR, LABP=>JLBP
      USE rdensity_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rintdens_I
      USE rinti_I
      USE rint_I
      USE rintdensvec_I
      USE cre_I
      USE getyn_I
      USE natorbnew_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*24, INTENT(IN) :: NAME
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(:,:,:), pointer :: DINT1VEC
      REAL(DOUBLE), DIMENSION(:,:), pointer   :: DENS1VEC
      REAL(DOUBLE), DIMENSION(:), pointer     :: DENSFIT
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: DINT1
      REAL(DOUBLE) :: AU2FM
      LOGICAL :: YES
      INTEGER :: I, J, L, ncount1
!-----------------------------------------------
!
      AU2FM = 52917.72083d0     ! JE ADD

!   Allocate memory

      ALLOCATE( DENS1VEC(NVEC,NNNP),DINT1VEC(NNNW,NNNW,NNNP) )    ! JE ADD
      ALLOCATE( DENSFIT(NNNP) )                                   ! JE ADD

      DINT1VEC(:,:,:) = 0.D0                                      ! JE ADD
      DENS1VEC(:,:) = 0.D0                                        ! JE ADD 

!
!   Allocate storage for local arrays
!
      CALL ALLOC (DENS1,NVEC,'DENS1','RDENSITY_CAL')

      CALL STARTTIME (ncount1, 'RDENSITY_CAL')

!
!   Initialise
!
      DO I = 1,NVEC
         DENS1(I) = 0.0D00
      ENDDO
 
!
!   Calculate all integrals needed for the volume shift calc.
!
!      IF (VSH) THEN
      DO 5 I = 1,NW
        DO 4 J = 1,NW
          IF (NAK(I).EQ.NAK(J)) THEN
            DINT1(I,J) = RINTDENS(I,J)
            CALL RINTDENSVEC(I,J,DINT1VEC)
          ELSE
            DINT1(I,J) = 0.0D00
            DINT1VEC(I,J,:) = 0.0D00
          ENDIF
    4   CONTINUE
    5 CONTINUE
!
!
     CALL NATORBNEW(NAME,DENS1VEC,DINT1VEC)                        
!
!
!
!   Printouts
!

      write(35,*) '     r [au]          D(r)=r^2*rho(r)       rho(r)'
!                 rho(r)-rho_fit(r) within nucleus'
      write(35,*)
      DO 16 I = 1,NVEC
	   WRITE (35,344)IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2)
           DO 17 L = 1,NNNP     ! JE ADD
              IF (DENS1VEC(I,L).GT.0.0) THEN                           ! JE ADD
                    WRITE (35,331) R(L),R(L)*R(L)*DENS1VEC(I,L),  &    ! JE ADD
                         DENS1VEC(I,L) ! JE ADD
	      ENDIF
 17        CONTINUE             ! JE ADD
 16   CONTINUE
!
!   Dealloc
!
      CALL DALLOC (DENS1,'DENS1','RDENSITY_CAL')

      CALL STOPTIME (ncount1, 'RDENSITY_CAL')
      RETURN
!
  331 FORMAT (1PD20.10,1PD20.10,1PD20.10)
  344 FORMAT (1X,I3,5X,2A4)
  
      RETURN
      END SUBROUTINE RDENSITY_CAL
