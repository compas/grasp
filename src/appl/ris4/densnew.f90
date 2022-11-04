!***********************************************************************
!                                                                      *
SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
!       SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,                      &
!                          DINT4,DINT5,DINT6,DINT7,                     &
!                          DINT1VEC,DENS1VEC)
!                                                                      *
!   IF angular coefficients must be calculated                         *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.
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
      USE debug_C
      USE decide_C
      USE DEF_C
      USE eigv_C
      USE foparm_C
      USE grid_C
      USE JLABL_C
      USE npar_C
      USE orb_C
      USE blk_C
      USE prnt_C
      USE TEILST_C
      USE BUFFER_C
      USE ris_C
      USE syma_C
      USE prnt_C,           ONLY : NVEC
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
                                                        DINT5, DINT6, &
                                                        DINT7
      INTEGER, INTENT(IN) :: DOIT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL, TSHELL_S
      REAL(DOUBLE), DIMENSION(:), pointer :: EMT1, EMT2, EMT3, EMT4,   &
                                             EMT5, EMT6
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5,     &
                      ELEMNT6, ELEMNT7, CONTRI1, CONTRI2, CONTRI3,     &
                      CONTRI4, CONTRI5, CONTRI6, CONTRI7
      LOGICAL :: VSH, NUCDE, SMSSH, YES
      CHARACTER :: CNUM*11, CK*2
      INTEGER, DIMENSION(NNNW) :: IA_S
      INTEGER :: KA, IOPAR, INCOR, IC, LCNUM, ITJPOC, IR, IA, IB, I, J
      Integer :: L, LOC, NCONTR, LAB, II
!-----------------------------------------------
!
! DOIT: IF DOIT=1 angular coefficients will be stored after creation
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
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
!   Matrix elements are diagonal in J
!
          IF (ITJPO(IR) .EQ. ITJPOC) THEN
!
!   Initialise the accumulator
!
            ELEMNT1 = 0.0D00
            ELEMNT2 = 0.0D00
            ELEMNT3 = 0.0D00
            ELEMNT4 = 0.0D00
            ELEMNT5 = 0.0D00
            ELEMNT6 = 0.0D00
            ELEMNT7 = 0.0D00
!
!   Call the MCT package to compute T coefficients
!
           CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
!GG            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)
            IF (IA .NE. 0) THEN
              IF (IA .EQ. IB) THEN
                    NCONTR = 0
                DO 8 IA = 1,NW
                  IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                    NCONTR = NCONTR + 1
                    TSHELL_S(NCONTR) = TSHELL(IA)
                    IA_S(NCONTR) = IA
                    ELEMNT1 = ELEMNT1 + DINT1(IA,IA)*TSHELL(IA)
                    ELEMNT2 = ELEMNT2 + DINT2(IA,IA)*TSHELL(IA)
                    ELEMNT3 = ELEMNT3 + DINT3(IA,IA)*TSHELL(IA)
                    ELEMNT4 = ELEMNT4 + DINT4(IA,IA)*TSHELL(IA)
                    ELEMNT5 = ELEMNT5 + DINT5(IA,IA)*TSHELL(IA)
                    ELEMNT6 = ELEMNT6 + DINT6(IA,IA)*TSHELL(IA)
                    ELEMNT7 = ELEMNT7 + DINT7(IA,IA)*TSHELL(IA)
                  ENDIF
    8           CONTINUE
                IF (DOIT.EQ.1) WRITE(50) IC,IR,NCONTR
                DO I = 1,NCONTR
                  LAB = IA_S(I)*(KEY + 1)
                  IF (DOIT.EQ.1) WRITE(50) TSHELL_S(I),LAB
                END DO
              ELSE
                IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                  IF (NAK(IA).EQ.NAK(IB)) THEN
                    IF (DOIT.EQ.1) WRITE(50) IC,IR,1
                    LAB = IA*KEY + IB
                    IF (DOIT.EQ.1) WRITE(50) TSHELL(1),LAB
                    ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL(1)
                    ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL(1)
                    ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL(1)
                    ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL(1)
                    ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL(1)
                    ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL(1)
                    ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL(1)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            DO 9 J = 1,NVEC
              LOC = (J-1)*NCF
              CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
              CONTRI2 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT2
              CONTRI3 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT3
              CONTRI4 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT4
              CONTRI5 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT5
              CONTRI6 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT6
              CONTRI7 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT7
              IF (IR.NE.IC) THEN
                CONTRI1 = 2.0D00 * CONTRI1
                CONTRI2 = 2.0D00 * CONTRI2
                CONTRI3 = 2.0D00 * CONTRI3
                CONTRI4 = 2.0D00 * CONTRI4
                CONTRI5 = 2.0D00 * CONTRI5
                CONTRI6 = 2.0D00 * CONTRI6
                CONTRI7 = 2.0D00 * CONTRI7
              ENDIF
              DENS1(J) = DENS1(J) + CONTRI1
              DENS2(J) = DENS2(J) + CONTRI2
              DENS3(J) = DENS3(J) + CONTRI3
              DENS4(J) = DENS4(J) + CONTRI4
              DENS5(J) = DENS5(J) + CONTRI5
              DENS6(J) = DENS6(J) + CONTRI6
              DENS7(J) = DENS7(J) + CONTRI7
    9       CONTINUE
          ENDIF
   12   CONTINUE
   13 CONTINUE

       IF (DOIT.EQ.1) WRITE(50) -1
!
! Empty the buffer and close file
      IF (DOIT.EQ.1) CLOSE(50)
!
!   Deallocate storage for the arrays in BUFFER
      CALL ALCBUF (3)
      RETURN
      END SUBROUTINE DENSNEW
