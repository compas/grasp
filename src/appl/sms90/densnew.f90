!***********************************************************************
!                                                                      *
      SUBROUTINE DENSNEW(DINT1, DINT2, DINT3, DINT4, DINT5, DINT6) 
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:42:57   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE debug_C
      USE decide_C
      USE DEF_C 
      USE eigv_C
      USE foparm_C
      USE grid_C
      USE JLABL_C 
      USE npar_C
      USE orb_C
      USE prnt_C 
      USE TEILST_C 
      USE BUFFER_C 
      USE SMS1_C 
      USE syma_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I 
      USE convrt_I 
      USE getyn_I
      USE itjpo_I 
      USE onescalar_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2, &
                                                        DINT3, DINT4, &
                                                        DINT5, DINT6
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KA, IOPAR, INCOR, IC, LCNUM, ITJPOC, IR, IA, IB, J, LOC
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL 
      REAL(DOUBLE), DIMENSION(:), pointer :: EMT1, EMT2, EMT3, EMT4, EMT5, EMT6
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5, ELEMNT6, &
         CONTRI1, CONTRI2, CONTRI3, CONTRI4, CONTRI5, CONTRI6 
      LOGICAL :: VSH, NUCDE, SMSSH, YES 
      CHARACTER :: CNUM*11, CK*2 
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      KA = 0 
      IOPAR = 1 
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
!   Matrix elements are diagonal in J
!
            IF (ITJPO(IR) /= ITJPOC) CYCLE  
!
!   Initialise the accumulator
!
            ELEMNT1 = 0.0D00 
            ELEMNT2 = 0.0D00 
            ELEMNT3 = 0.0D00 
            ELEMNT4 = 0.0D00 
            ELEMNT5 = 0.0D00 
            ELEMNT6 = 0.0D00 
!
!   Call the MCT package to compute T coefficients
!
           CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
!            CALL TNSRJJ (KA, IOPAR, IC, IR, IA, IB, TSHELL) 
 
            IF (IA /= 0) THEN 
               IF (IA == IB) THEN 
                  DO IA = 1, NW 
                     IF (ABS(TSHELL(IA)) <= CUTOFF) CYCLE  
                     ELEMNT1 = ELEMNT1 + DINT1(IA,IA)*TSHELL(IA) 
                     ELEMNT2 = ELEMNT2 + DINT2(IA,IA)*TSHELL(IA) 
                     ELEMNT3 = ELEMNT3 + DINT3(IA,IA)*TSHELL(IA) 
                     ELEMNT4 = ELEMNT4 + DINT4(IA,IA)*TSHELL(IA) 
                     ELEMNT5 = ELEMNT5 + DINT5(IA,IA)*TSHELL(IA) 
                     ELEMNT6 = ELEMNT6 + DINT6(IA,IA)*TSHELL(IA) 
                  END DO 
               ELSE 
                  IF (ABS(TSHELL(1)) > CUTOFF) THEN 
                     IF (NAK(IA) == NAK(IB)) THEN 
                        ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL(1) 
                        ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL(1) 
                        ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL(1) 
                        ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL(1) 
                        ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL(1) 
                        ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL(1) 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
            DO J = 1, NVEC 
               LOC = (J - 1)*NCF 
               CONTRI1 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT1 
               CONTRI2 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT2 
               CONTRI3 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT3 
               CONTRI4 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT4 
               CONTRI5 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT5 
               CONTRI6 = EVEC(IC + LOC)*EVEC(IR + LOC)*ELEMNT6 
               IF (IR /= IC) THEN 
                  CONTRI1 = 2.0D00*CONTRI1 
                  CONTRI2 = 2.0D00*CONTRI2 
                  CONTRI3 = 2.0D00*CONTRI3 
                  CONTRI4 = 2.0D00*CONTRI4 
                  CONTRI5 = 2.0D00*CONTRI5 
                  CONTRI6 = 2.0D00*CONTRI6 
               ENDIF 
               DENS1(J) = DENS1(J) + CONTRI1 
               DENS2(J) = DENS2(J) + CONTRI2 
               DENS3(J) = DENS3(J) + CONTRI3 
               DENS4(J) = DENS4(J) + CONTRI4 
               DENS5(J) = DENS5(J) + CONTRI5 
               DENS6(J) = DENS6(J) + CONTRI6 
            END DO 
         END DO 
      END DO 
!
!   Deallocate storage for the arrays in BUFFER
!
      CALL ALCBUF (3) 
      RETURN  
      END SUBROUTINE DENSNEW 
