!***********************************************************************
!                                                                      *
      SUBROUTINE SMS
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the  sms parameter, the electron density at the    *
!   origin.
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, GETYN          *
!                        ITJPO, RKCO, TNSRJJ                           *
!               [SMS92]: RINTISO, RINTDENS, VINTI                      *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!***********************************************************************
!...Created by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
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
      USE sms1_C
      USE jlabl_C, LABJ=> JLBR, LABP=>JLBP
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rintdens_I
      USE rinti_I
      USE rint_I
      USE vinti_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW)      :: TSHELL
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: VINT
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: DINT1, DINT2, DINT3, &
                                            DINT4, DINT5, DINT6
      CHARACTER*11 :: CNUM
      CHARACTER*2  :: CK
      LOGICAL ::  GETYN,VSH,NUCDE,SMSSH,YES,AVAIL
      INTEGER :: I, J
!-----------------------------------------------
!
!   Allocate storage for local arrays
!
      CALL ALLOC (SMSC, NVEC,'SMSC', 'SMS')
      CALL ALLOC (DENS1,NVEC,'DENS1','SMS')
      CALL ALLOC (DENS2,NVEC,'DENS2','SMS')
      CALL ALLOC (DENS3,NVEC,'DENS3','SMS')
      CALL ALLOC (DENS4,NVEC,'DENS4','SMS')
      CALL ALLOC (DENS5,NVEC,'DENS5','SMS')
      CALL ALLOC (DENS6,NVEC,'DENS6','SMS')
!
!   Initialise
!
      DO I = 1,NVEC
         SMSC(I)  = 0.0D00
         DENS1(I) = 0.0D00
         DENS2(I) = 0.0D00
         DENS3(I) = 0.0D00
         DENS4(I) = 0.0D00
         DENS5(I) = 0.0D00
         DENS6(I) = 0.0D00
      ENDDO

      VSH   = .TRUE.
      SMSSH = .TRUE.
!
!   Calculate all integrals needed for the volume shift calc.
!
      IF (VSH) THEN
         DO 5 I = 1,NW
            DO 4 J = 1,NW
               IF (NAK(I).EQ.NAK(J)) THEN
                  DINT1(I,J) = RINTDENS(I,J)
                  DINT2(I,J) = RINTI(I,J,1)
                  DINT3(I,J) = RINT(I,J,1)
                  DINT4(I,J) = RINT(I,J,2)
                  DINT5(I,J) = RINT(I,J,-1)
                  DINT6(I,J) = RINT(I,J,-2)
               ELSE
                  DINT1(I,J) = 0.0D00
                  DINT2(I,J) = 0.0D00
                  DINT3(I,J) = 0.0D00
                  DINT4(I,J) = 0.0D00
                  DINT5(I,J) = 0.0D00
                  DINT6(I,J) = 0.0D00
               ENDIF
    4       CONTINUE
    5    CONTINUE
      ENDIF
!
!   Calculate and save the Vinti integrals
!
      IF (SMSSH) THEN
         DO 7 I = 1,NW
            DO 6 J = 1,NW
               IF (I.NE.J) THEN
                  VINT(I,J) = VINTI(I,J)
               ELSE
                  VINT(I,J) = 0.0D00
               ENDIF
    6       CONTINUE
    7    CONTINUE
      ENDIF
!
!   See if the appropriate angular data is available. If so,
!   then read the angular files and perform the calculation.
!
      IF (SMSSH) THEN
         CALL SETMCP(AVAIL)
         IF (AVAIL) THEN
            CALL SMSMCP(VINT)
         ELSE
            CALL SMSNEW(VINT)
         ENDIF
      ENDIF
      IF (VSH) THEN
         CALL SETMCP(AVAIL)
         IF (AVAIL) THEN
            CALL DENSMCP(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)
         ELSE
            CALL DENSNEW(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6)
         ENDIF
      ENDIF
!
!   Printouts
!
      WRITE (24,301) CUTOFF
      WRITE (24,302)
      DO 14 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        SMSC(I)
   14 CONTINUE
      WRITE (24,307)
      DO 16 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS1(I)
   16 CONTINUE
      WRITE (24,308)
      DO 18 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS2(I)
   18 CONTINUE
      WRITE (24,309)
      DO 20 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS3(I)
   20 CONTINUE
      WRITE (24,310)
      DO 22 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS4(I)
   22 CONTINUE
      WRITE (24,311)
      DO 24 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS5(I)
   24 CONTINUE
      WRITE (24,312)
      DO 26 I = 1,NVEC
         WRITE (24,303) IVEC(I),LABJ(IATJPO(I)),LABP((IASPAR(I)+3)/2), &
                        DENS6(I)
   26 CONTINUE
!
!   Dealloc
!
      CALL DALLOC (SMSC, 'SMSC', 'SMS')
      CALL DALLOC (DENS1,'DENS1','SMS')
      CALL DALLOC (DENS2,'DENS2','SMS')
      CALL DALLOC (DENS3,'DENS3','SMS')
      CALL DALLOC (DENS4,'DENS4','SMS')
      CALL DALLOC (DENS5,'DENS5','SMS')
      CALL DALLOC (DENS6,'DENS6','SMS')
!
  301 FORMAT (//' CUTOFF set to ',1PD22.15)
  302 FORMAT (//' Level  J Parity  Specific mass shift (au) '/)
  303 FORMAT (1X,I3,5X,2A4,3X,D20.10)
  307 FORMAT (//' Electron density in atomic units'                    &
              //' Level  J Parity',8X,'DENS (a.u.)'/)
  308 FORMAT (//' Kinetic energy '                                     &
              //' Level  J Parity',8X,'T (a.u.)'/)
  309 FORMAT (//' Radial expectationvalue'                             &
              //' Level  J Parity',8X,'<r> (a.u.)'/)
  310 FORMAT (//' Radial expectationvalue'                             &
              //' Level  J Parity',8X,'<r2> (a.u.)'/)
  311 FORMAT (//' Radial expectationvalue'                             &
              //' Level  J Parity',8X,'<r-1> (a.u.)'/)
  312 FORMAT (//' Radial expectationvalue'                             &
              //' Level  J Parity',8X,'<r-2> (a.u.)'/)

!
      RETURN
      END SUBROUTINE SMS
