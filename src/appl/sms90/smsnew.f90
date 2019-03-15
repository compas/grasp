!***********************************************************************
!                                                                      *
      SUBROUTINE SMSNEW(VINT)
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
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                         Last revision:    Nov 2017   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:48:15   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE BUFFER_C
      USE TEILST_C
      USE SMS1_C
      USE debug_C
      USE decide_C
      USE def_C
      USE eigv_C
      USE foparm_C
      USE grid_C
      USE npar_C
      USE prnt_C
      USE syma_C
      USE orb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cord_I
      USE alcbuf_I
      USE convrt_I
      USE itjpo_I
      USE rkco_gg_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: VINT(NNNW,NNNW)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: INCOR, IC, LCNUM, ITJPOC, IR, I, IIA, IIB, IIC, IID, K, J, LOC
      REAL(DOUBLE) :: VCOEFF, CONTRI
      LOGICAL :: GETYN, VSH, NUCDE, SMSSH, YES
      CHARACTER :: CNUM*11, CK*2
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
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
      DO IC = 1, NCF
!
!   Output IC on the screen to show how far the calculation has preceede
!
         CALL CONVRT (IC, CNUM, LCNUM)
         IF (MOD(IC,10) == 0) WRITE (6, *) 'Column '//CNUM(1:LCNUM)//&
            ' complete;'
!
         ITJPOC = ITJPO(IC)
         DO IR = IC, NCF
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs
!
!   Initialize
!
!
!   Matrix elements are diagonal in J
!
            IF (ITJPO(IR) /= ITJPOC) CYCLE
            NVCOEF = 0
!GG            CALL RKCO (IC, IR, COR, CORD, INCOR)
            CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
!
            DO I = 1, NVCOEF
               VCOEFF = COEFF(I)
               IF (ABS(VCOEFF) <= CUTOFF) CYCLE
               IIA = LABEL(1,I)
               IIB = LABEL(2,I)
               IIC = LABEL(3,I)
               IID = LABEL(4,I)
               K = LABEL(5,I)
!
!   Only K = 1   LABEL(5,I) .EQ. 1
!
               IF (LABEL(5,I) /= 1) CYCLE
               IF (LDBPA(2)) WRITE (99, 309) K, IC, IR, NP(IIA), NH(IIA), NP(&
                  IIB), NH(IIB), NP(IIC), NH(IIC), NP(IID), NH(IID), VCOEFF
               DO J = 1, NVEC
                  LOC = (J - 1)*NCF
                  CONTRI = -EVEC(IC + LOC)*EVEC(IR + LOC)*VCOEFF*VINT(LABEL(1,I&
                     ),LABEL(3,I))*VINT(LABEL(2,I),LABEL(4,I))
                  IF (IR /= IC) CONTRI = 2.0D00*CONTRI
                  SMSC(J) = SMSC(J) + CONTRI
               END DO
            END DO
         END DO
      END DO
!
!   Deallocate storage for the arrays in BUFFER
!
      CALL ALCBUF (3)
      RETURN
  309 FORMAT(' V^[(',1I2,')]_[',1I3,',',1I3,']',' (',1I2,1A2,',',1I2,1A2,';',1I&
         2,1A2,',',1I2,1A2,') = ',1P,D19.12)
      RETURN
!
      END SUBROUTINE SMSNEW
