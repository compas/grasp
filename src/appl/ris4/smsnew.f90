!***********************************************************************
!                                                                      *
      SUBROUTINE SMSNEW(DOIT,VINT,VINT2)
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the  sms parameter, the electron density at the    *
!   origin.
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
!                        ITJPO,                                        *
!               [SMS92]: RINTISO, RINTDENS, VINTI                      *
!               [LIBRANG]: RKCO_GG                                     *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Mai. 2012                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE debug_C
      USE prnt_C
      USE orb_C
      USE BUFFER_C
      USE eigv_C
      USE ris_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
      INTEGER, INTENT(IN) :: DOIT
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL     :: CORD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*11 :: CNUM
      CHARACTER*4  :: JLBL,LABJ,LABP
      CHARACTER*2  :: CK
      LOGICAL      :: GETYN,FIRSTT,VSH,NUCDE,SMSSH,YES
      LOGICAL      :: LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      REAL(DOUBLE) :: CONTRI, CONTRIK1, VCOEFF
      INTEGER      :: I,J,K,LAB,LOC,IC,IR,ITJPOC,IIA,IIB,IIC,IID,INCOR
      INTEGER      :: LCNUM
!-----------------------------------------------
!
      INCOR = 1
!
!   Allocate storage for the arrays in BUFFER
!
      CALL ALCBUF (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
      DO 13 IC = 1,NCF
!
!   Output IC on the screen to show how far the calculation has preceede
!
        CALL CONVRT (IC,CNUM,LCNUM)
        if (mod(IC,100).eq.0) then
          PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
        end if
!
        ITJPOC = ITJPO (IC)
        DO 12 IR = IC,NCF
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs
!
!   Initialize
!
!
!   Matrix elements are diagonal in J
!
          IF (ITJPO(IR) .EQ. ITJPOC) THEN
            NVCOEF = 0
!GG            CALL RKCO (IC,IR,COR,CORD,INCOR)
            CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
!
            DO 11 I = 1,NVCOEF
              VCOEFF = COEFF(I)
              IF (ABS (VCOEFF) .GT. CUTOFF) THEN
                IIA = LABEL(1,I)
                IIB = LABEL(2,I)
                IIC = LABEL(3,I)
                IID = LABEL(4,I)
                K   = LABEL(5,I)
!
!   Only K = 1   LABEL(5,I) .EQ. 1
!
                IF (LABEL(5,I) .EQ. 1) THEN
                  IF (LDBPA(2)) THEN
                    WRITE (99,309) K,IC,IR,                            &
                    NP(IIA),NH(IIA),NP(IIB),NH(IIB),                   &
                    NP(IIC),NH(IIC),NP(IID),NH(IID),VCOEFF
                  ENDIF
                  IF(DOIT.EQ.1) WRITE(51) IC,IR
!** Storage sequence
                  LAB  = ((IIA*KEY + IIC)*KEY+IIB)*KEY+IID
                  IF(DOIT.EQ.1) THEN
                    WRITE(51) VCOEFF,LAB
                  ENDIF
!***
                  DO 10 J = 1,NVEC
                    LOC = (J-1)*NCF
                    CONTRIK1 = - EVEC(IC+LOC)*EVEC(IR+LOC)             &
                           * VCOEFF                                    &
                           * VINT (LABEL(1,I),LABEL(3,I))              &
                           * VINT (LABEL(2,I),LABEL(4,I))
                    CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)               &
                           * VCOEFF                                    &
                           * ( VINT2(LABEL(1,I),LABEL(3,I))            &
                           * VINT (LABEL(2,I),LABEL(4,I))              &
                           + VINT2(LABEL(2,I),LABEL(4,I))              &
                           * VINT (LABEL(1,I),LABEL(3,I)) )/2.0D00
                    IF (IR.NE.IC) THEN
                      CONTRI = 2.0D00 * CONTRI
                      CONTRIK1 = 2.0D00 * CONTRIK1
                    ENDIF
                    SMSC1(J) = SMSC1(J) + CONTRIK1
                    SMSC2(J) = SMSC2(J) + CONTRI
   10             CONTINUE
                ENDIF
              ENDIF
   11       CONTINUE
          ENDIF
   12   CONTINUE
   13 CONTINUE
      IF(DOIT.EQ.1) WRITE(51) -1
!
! Empty the buffer and close file
      IF(DOIT.EQ.1) CLOSE(51)
!
!   Deallocate storage for the arrays in BUFFER
!
      CALL ALCBUF (3)
      RETURN
  309 FORMAT (' V^[(',1I2,')]_[',1I3,',',1I3,']',                      &
         ' (',1I2,1A2,',',1I2,1A2,';',                                 &
              1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
!
      END SUBROUTINE SMSNEW
