!***********************************************************************
!                                                                      *
      SUBROUTINE RIS_CAL (NAME)
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
      USE ris_C
      USE jlabl_C, LABJ=> JLBR, LABP=>JLBP
      USE ris_C
      USE eigv_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rintdens_I
      USE rinti_I
      USE rint_I
      USE rintdensvec_I
      USE rinti_nms_I
      USE cre_I
      USE vinti_I
      USE rint_sms2_I
      USE rint_sms3_I
      USE angdata_I
      USE getyn_I
      USE densread_I
      USE densnew_I
      USE smsread_I
      USE smsnew_I
      USE edensityfit_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*24, INTENT(IN) :: NAME
!      REAL(DOUBLE), INTENT(IN) :: DR2
!      INTEGER, INTENT(IN)      :: NOPAR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(:,:,:), pointer :: DINT1VEC
      REAL(DOUBLE), DIMENSION(:,:), pointer   :: DENS1VEC
      REAL(DOUBLE), DIMENSION(:), pointer     :: DENSFIT
      REAL(DOUBLE), DIMENSION(:,:), pointer     :: FMAT
      REAL(DOUBLE), DIMENSION(:), pointer     :: RHO
      REAL(DOUBLE), DIMENSION(:), pointer     :: RES


      REAL(DOUBLE), DIMENSION(NNNW)      :: TSHELL
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: VINT, VINT2
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: DINT1, DINT2, DINT3, &
                                            DINT4, DINT5, DINT6, &
                                            DINT7
      REAL(DOUBLE) :: AU2FM, RCRE, HzSMSu, HzNMSu, FO90
      CHARACTER*11 :: CNUM
      CHARACTER*2  :: CK
      LOGICAL :: VSH, NUCDE, SMSSH, YES, YES2, AVAIL_TB, AVAIL_OB
      INTEGER :: I, J, L, ncount1, DOIT_OB, DOIT_TB, NRNUC
!-----------------------------------------------
!

!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!      PARAMETER (CUTOFF = 1.0D-10)

      write(*,*) '-------------------------------'
      write(*,*) 'RIS_CAL: Execution Begins ...'
      write(*,*) '-------------------------------'

!     Determine point end point NRNUC for density fit
      FO90 = PARM(1)+2.d0*LOG(3.d0)*PARM(2)
      !PARM(1) - Fermi parameter c
      !PARM(2) - Fermi parameter a
      DO I=1,NNNP
         IF(R(I+1).GT.FO90.AND.R(I).LE.FO90) THEN
            NRNUC = I+1
            exit
         END IF
      END DO


!   Allocate memory

      ALLOCATE( DENS1VEC(NVEC,NRNUC),DINT1VEC(NNNW,NNNW,NRNUC) )    ! JE ADD
      ALLOCATE( DENSFIT(NRNUC) )                                   ! JE ADD

      DINT1VEC(:,:,:) = 0.D0                                      ! JE ADD
      DENS1VEC(:,:) = 0.D0                                        ! JE ADD

      ALLOCATE( FMAT(6,NVEC),RHO(NVEC),RES(NVEC) )                ! JE ADD

!
!   Allocate storage for local arrays
!
      CALL ALLOC (SMSC1, NVEC,'SMSC1', 'RIS_CAL')
      CALL ALLOC (SMSC2, NVEC,'SMSC2', 'RIS_CAL')
      CALL ALLOC (DENS1,NVEC,'DENS1','RIS_CAL')
      CALL ALLOC (DENS2,NVEC,'DENS2','RIS_CAL')
      CALL ALLOC (DENS3,NVEC,'DENS3','RIS_CAL')
      CALL ALLOC (DENS4,NVEC,'DENS4','RIS_CAL')
      CALL ALLOC (DENS5,NVEC,'DENS5','RIS_CAL')
      CALL ALLOC (DENS6,NVEC,'DENS6','RIS_CAL')
      CALL ALLOC (DENS7,NVEC,'DENS7','RIS_CAL')

      CALL STARTTIME (ncount1, 'RIS_CAL')

!
!   Initialise
!
      DO I = 1,NVEC
         SMSC1(I)  = 0.0D00
         SMSC2(I)  = 0.0D00
         DENS1(I) = 0.0D00
         DENS2(I) = 0.0D00
         DENS3(I) = 0.0D00
         DENS4(I) = 0.0D00
         DENS5(I) = 0.0D00
         DENS6(I) = 0.0D00
         DENS7(I) = 0.0D00
      ENDDO

!      VSH   = .TRUE.
!      SMSSH = .TRUE.
!
!   Calculate all integrals needed for the volume shift calc.
!
      !      IF (VSH) THEN

      write(*,*) ' Compute higher order field shift electronic factors?'
      YES2 = GETYN ()

      DO 5 I = 1,NW
        DO 4 J = 1,NW
          IF (NAK(I).EQ.NAK(J)) THEN
            DINT1(I,J) = RINTDENS(I,J)
            !GG                 DINT2(I,J) = RINTI(I,J,1)
            IF(YES2) THEN
               CALL RINTDENSVEC(I,J,DINT1VEC,NRNUC)
            END IF
            CALL RINTI_NMS(I,J,DINT2(I,J),DINT7(I,J))
            DINT3(I,J) = RINT(I,J,1)
            DINT4(I,J) = RINT(I,J,2)
            DINT5(I,J) = RINT(I,J,-1)
            DINT6(I,J) = RINT(I,J,-2)
          ELSE
            DINT1(I,J) = 0.0D00
            DINT1VEC(I,J,:) = 0.0D00
            DINT2(I,J) = 0.0D00
            DINT3(I,J) = 0.0D00
            DINT4(I,J) = 0.0D00
            DINT5(I,J) = 0.0D00
            DINT6(I,J) = 0.0D00
            DINT7(I,J) = 0.0D00
          ENDIF
    4   CONTINUE
    5 CONTINUE
!      ENDIF
!
!   Calculate and save the Vinti integrals
!
!      IF (SMSSH) THEN
      DO 7 I = 1,NW
        DO 6 J = 1,NW
          IF (I.NE.J) THEN
            RCRE = CRE(NAK(I),1,NAK(J))
            IF (DABS(RCRE) .GT. CUTOFF) THEN
              VINT  (I,J) = VINTI(I,J)
              VINT2(I,J) = VINT(I,J)                                   &
                         + RINT_SMS2(I,J)/RCRE                         &
                         + RINT_SMS3(I,J)
            ELSE
             VINT (I,J) = 0.0D00
             VINT2(I,J) = 0.0D00
            ENDIF
          ELSE
            VINT (I,J) = 0.0D00
            VINT2(I,J) = 0.0D00
          ENDIF
    6   CONTINUE
    7 CONTINUE
!      ENDIF
!
!   See if the appropriate angular data is available. If so,
!   then read the angular files and perform the calculation.
!   If not, the user can choose to save them or not
!
      DOIT_OB = 0
      DOIT_TB = 0
      CALL ANGDATA(NAME,AVAIL_OB,1)
      CALL ANGDATA(NAME,AVAIL_TB,2)
      IF ((.NOT. AVAIL_OB) .AND. (.NOT. AVAIL_TB)) THEN
        PRINT *,' Save ang. coefficients of one- and two-body op.?'
        YES = GETYN ()
        PRINT *
!C        IF(.NOT. YES) THEN
!C          PRINT *,' Save ang. coefficients of one-body op. ONLY?'
!C          YES = GETYN ()
!C          PRINT *
!C          IF(.NOT. YES) THEN
!C            PRINT *,' Save ang. coefficients of two-body op. ONLY?'
!C            YES = GETYN ()
!C            PRINT *
!C            IF(YES) DOIT_TB = 1
!C          ELSE
!C            DOIT_OB = 1
!C          ENDIF
!C        ELSE
       IF(YES) THEN
          DOIT_OB = 1
          DOIT_TB = 1
        ENDIF
      ELSEIF (.NOT. AVAIL_OB) THEN
        PRINT *,' Save ang. coefficients of one-body op. ?'
        YES = GETYN ()
        PRINT *
        IF(YES) DOIT_OB = 1
      ELSEIF (.NOT. AVAIL_TB) THEN
        PRINT *,' Save ang. coefficients of two-body op. ?'
        YES = GETYN ()
        PRINT *
        IF(YES) DOIT_TB = 1
      ENDIF
      IF (AVAIL_OB) THEN
        J = INDEX(NAME,' ')
        OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN'       &
            ,FORM='UNFORMATTED')
!    Read angular data from file and compute matrix elements
        IF(YES2) THEN
           CALL DENSREAD_SELTZ(DINT1,DINT2,DINT3, & ! JE ADD
                DINT4,DINT5,DINT6, & ! JE ADD
                DINT7,DINT1VEC,DENS1VEC,NRNUC) ! JE ADD
        ELSE
           CALL DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
        END IF

        !        CALL DENSREAD(DINT1,DINT2,DINT3,             &   ! JE ADD
        !                      DINT4,DINT5,DINT6,             &   ! JE ADD
        !                      DINT7,DINT1VEC,DENS1VEC)           ! JE ADD
     ELSE
!Check if the user wants to save one-body ang. coeff.
        IF(DOIT_OB .EQ.1) THEN
          J = INDEX(NAME,' ')
          OPEN(UNIT=50,FILE = NAME(1:J-1)//'.IOB',STATUS='UNKNOWN'     &
            ,FORM='UNFORMATTED')
        ENDIF
        IF(YES2) THEN
           CALL DENSNEW_SELTZ(DOIT_OB,DINT1,DINT2,DINT3, & ! JE ADD
                DINT4,DINT5,DINT6,DINT7, & ! JE ADD
                DINT1VEC,DENS1VEC,NRNUC) ! JE ADD
        ELSE
           CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,DINT4, &
                DINT5,DINT6,DINT7)
        END IF

!        CALL NATORBNEW(NAME)                            ! JE ADD
!        CALL DENSNEW(DOIT_OB,DINT1,DINT2,DINT3,     &   ! JE ADD
!             DINT4,DINT5,DINT6,DINT7,               &   ! JE ADD
!             DINT1VEC,DENS1VEC)                         ! JE ADD
      ENDIF
!
!
!
      IF (AVAIL_TB) THEN
        J = INDEX(NAME,' ')
        OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN'       &
          ,FORM='UNFORMATTED')
        CALL SMSREAD(VINT,VINT2)
      ELSE
! Check if the user wants to save two-body ang. coeff.
        IF(DOIT_TB.EQ.1) THEN
          J = INDEX(NAME,' ')
          OPEN(UNIT=51,FILE = NAME(1:J-1)//'.ITB',STATUS='UNKNOWN'     &
             ,FORM='UNFORMATTED')
        ENDIF
         CALL SMSNEW(DOIT_TB,VINT,VINT2)
      ENDIF
!
!   Printouts
      !

      WRITE (24,'(a24,i3)') 'Number of eigenvalues: ',NVEC
      WRITE(24,*)
      WRITE(24,*)
      WRITE (24,309)
      DO I = 1, NVEC
         WRITE (24,323) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              EAV+EVAL(I)
      END DO
      WRITE (24,308)
      DO 18 I = 1,NVEC
         WRITE(24,314)
         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              DENS7(I),(DENS2(I)-DENS7(I)),DENS2(I)
         HzNMSu=DENS7(I)*AUMAMU*AUCM*CCMS
!         if (dabs(HzNMSu) .LE. (10**8)) THEN
!             WRITE (24,315) HzNMSu*(1.0D-06),
!     :      (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-06),
!     :       DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
!         else
         WRITE (24,316) HzNMSu*(1.0D-09), &
              (DENS2(I)-DENS7(I))*AUMAMU*AUCM*CCMS*(1.0D-09), &
              DENS2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
!         endif
 18   CONTINUE

      WRITE (24,302)
      DO 14 I = 1,NVEC
         WRITE(24,314)
         WRITE (24,313) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
              SMSC1(I),SMSC2(I)-SMSC1(I),SMSC2(I)
         HzSMSu=SMSC1(I)*AUMAMU*AUCM*CCMS
         if (dabs(HzSMSu) .LE. (10**8)) THEN
            WRITE (24,315) HzSMSu*(1.0D-06), &
                 (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-06), &
                 SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-06)
         else
            WRITE (24,316) HzSMSu*(1.0D-09), &
                 (SMSC2(I)-SMSC1(I))*AUMAMU*AUCM*CCMS*(1.0D-09), &
                 SMSC2(I)*AUMAMU*AUCM*CCMS*(1.0D-09)
         endif
 14   CONTINUE
      IF(YES2) THEN
         WRITE (24,307)
         DO 16 I = 1,NVEC
            CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT(:,I),RHO(I),RES(I))  ! JE ADD
            WRITE (24,303) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),RHO(I)
 16      CONTINUE
         WRITE (24,317)
         DO 17 I = 1,NVEC
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            WRITE (24,324) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),FMAT(1,I),FMAT(2,I), &
                 FMAT(3,I),FMAT(4,I),RES(I)
 17      CONTINUE
         WRITE (24,327)
         DO 19 I = 1,NVEC
            !CALL EDENSITYFIT(R,DENS1VEC(I,:),Z,PARM,NRNUC,FMAT,RHO,RES)  ! JE ADD
            WRITE (24,334) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),FMAT(5,I),FMAT(6,I)
 19      CONTINUE
         write(24,*)
      ELSE
         WRITE (24,306)
         DO 36 I = 1,NVEC
            WRITE (24,303) IVEC(I),LABJ(IATJPO(I)), &
                 LABP((IASPAR(I)+3)/2),DENS1(I)
 36      CONTINUE
      END IF
!     WRITE (24,309)
!      DO 20 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS3(I)
!     20 CONTINUE
!     WRITE (24,310)
!     DO 22 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS4(I)
!     22 CONTINUE
!     WRITE (24,311)
!     DO 24 I = 1,NVEC
!         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS5(I)
!     24 CONTINUE
!     WRITE (24,312)
!     DO 26 I = 1,NVEC
!     WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2),
!     :                  DENS6(I)
!     26 CONTINUE
!
!   Dealloc
!
      CALL DALLOC (SMSC1,'SMSC1','RIS_CAL')
      CALL DALLOC (SMSC2,'SMSC2','RIS_CAL')
      CALL DALLOC (DENS1,'DENS1','RIS_CAL')
      CALL DALLOC (DENS2,'DENS2','RIS_CAL')
      CALL DALLOC (DENS3,'DENS3','RIS_CAL')
      CALL DALLOC (DENS4,'DENS4','RIS_CAL')
      CALL DALLOC (DENS5,'DENS5','RIS_CAL')
      CALL DALLOC (DENS6,'DENS6','RIS_CAL')
      CALL DALLOC (DENS7,'DENS7','RIS_CAL')

      write(*,*) '-------------------------------'
      write(*,*) 'RIS_CAL: Execution Finished ...'
      write(*,*) '-------------------------------'

      CALL STOPTIME (ncount1, 'RIS_CAL')
      RETURN
!
301   FORMAT (//' CUTOFF set to ',1PD22.15)
302   FORMAT (//' Level  J Parity  Specific mass shift parameter')

!     Electron density values
303   FORMAT (1X,I3,5X,2A4,3X,D20.10)
304   FORMAT (1X,I3,5X,2A4,3X,2D20.10)

!     Electronic factor values
324   FORMAT (1X,I3,5X,2A4,3X,4D20.10,F10.4)

!     Field shift values
334   FORMAT (1X,I3,5X,2A4,3X,2D20.10)

      !     Electron density header
306   FORMAT (//' Electron density in atomic units' &
           //' Level  J Parity',8X,'DENS (a.u.)'/)

307   FORMAT (//' Level  J Parity  Electron density in atomic units' &
           //24X,'Dens. (a.u.)')

!     Electronic factors header
317   FORMAT (//' Level  J Parity  Field shift electronic factors and av&
           erage point discrepancy in fit' &
           //24X, &
           'F0 (GHz/fm^2)',6X,' F2 (GHz/fm^4)', &
           6X,' F4 (GHz/fm^6)',6X,' F6 (GHz/fm^8)' &
           6X,' Disc. (per mille)')

!     Electronic field shift header
327   FORMAT (//' Level  J Parity  Field shift electronic factors (corre&
           cted for varying density inside nucleus)' &
           //24X, &
           'F0VED0 (GHz/fm^2)   F0VED1 (GHz/fm^4)')

308   FORMAT (//' Level  J Parity  Normal mass shift parameter')

309   FORMAT (' Level  J Parity  Energy')

!     309 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r> (a.u.)'/)
!     310 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r2> (a.u.)'/)
!     311 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-1> (a.u.)'/)
!     312 FORMAT (//' Radial expectationvalue'
!     :        //' Level  J Parity',8X,'<r-2> (a.u.)'/)
313   FORMAT (1X,I3,5X,2A4,3x,3D20.10,'  (a.u.)')
333   FORMAT (1X,I3,5X,2A4,3x,3D20.10,'  (GHz u)')
323   FORMAT (1X,I3,5X,2A4,3x,D20.10,'  (a.u.)')
314   FORMAT (/,29X,'<K^1>',13X,'<K^2+K^3>',9X,'<K^1+K^2+K^3>')
315   FORMAT (20X,3D20.10 ,2X,'(MHz u)')
316   FORMAT (20X,3D20.10,2X,'(GHz u)')
      !
      RETURN
    END SUBROUTINE RIS_CAL
