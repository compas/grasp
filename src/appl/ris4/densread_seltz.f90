!***********************************************************************
!                                                                      *
!JE   SUBROUTINE DENSREAD(DINT1,DINT2,DINT3,DINT4,DINT5,DINT6,DINT7)
      SUBROUTINE DENSREAD_SELTZ(DINT1,DINT2,DINT3,                     &
                          DINT4,DINT5,DINT6,                           &
                          DINT7,DINT1VEC,DENS1VEC,NRNUC)
!                                                                      *
!   IF angular coefficients already exist                              *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.                              *
!                                                                      *
!   Written by C\'edric Naz\'e                                         *
!                                                                      *
!                                         Last revision: Feb. 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW, NNNP
      USE prnt_C
      USE ris_C
      USE orb_C
      USE eigv_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2,  &
                                          DINT3, DINT4, DINT5, DINT6,  &
                                          DINT7
      REAL(DOUBLE), DIMENSION(NVEC,NRNUC), INTENT(OUT)     :: DENS1VEC !  JE ADD
      REAL(DOUBLE), DIMENSION(NNNW,NNNW,NRNUC), INTENT(IN) :: DINT1VEC !  JE ADD
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL_R
      REAL(DOUBLE), DIMENSION(NRNUC) :: CONTRI1VEC, ELEMNT1VEC
      REAL(DOUBLE) :: ELEMNT1, ELEMNT2, ELEMNT3, ELEMNT4, ELEMNT5
      REAL(DOUBLE) :: ELEMNT6, ELEMNT7
      REAL(DOUBLE) :: CONTRI1, CONTRI2, CONTRI3, CONTRI4, CONTRI5
      REAL(DOUBLE) :: CONTRI6, CONTRI7
      INTEGER :: ICOLD, IROLD, IOS, IA, IB, IC, IR, I, J, L, LOC, LAB
      INTEGER :: NCOUNT, NRNUC
!-----------------------------------------------
!
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
!
      DENS1VEC(:,:) = 0.0D00                                     ! JE ADD

      ICOLD = 0
      IROLD = 0

      REWIND (50)
   16 READ (50,IOSTAT = IOS) IC,IR,NCOUNT

      IF (IOS .EQ. 0) THEN
!*      print*, 'ic', IC,IR,NCOUNT
!
!   Read successful; decode the labels of I(ab)
!
!      Initialise
        IF ((IC.NE.ICOLD).OR.(IR.NE.IROLD)) THEN
!*            ISPARC = ISPAR (IC)
!*            ITJPOC = ITJPO (IC)
!*            ITJPOR = ITJPO (IR)
!*            IDIFF = ITJPOC - ITJPOR
            ICOLD = IC
            IROLD = IR
!*        DO I = 1,NELMNT
          ELEMNT1 = 0.0D00
          ELEMNT2 = 0.0D00
          ELEMNT3 = 0.0D00
          ELEMNT4 = 0.0D00
          ELEMNT5 = 0.0D00
          ELEMNT6 = 0.0D00
          ELEMNT7 = 0.0D00
          ELEMNT1VEC(:) = 0.0D00                               ! JE ADD
!*        ENDDO
        ENDIF
        DO I= 1,NCOUNT
           READ(50) TSHELL_R(I),LAB
           IA = LAB/KEY
           IB = MOD(LAB,KEY)
!*           LOC = INDEX(I)
           ELEMNT1 = ELEMNT1 + DINT1(IA,IB)*TSHELL_R(I)
           DO 21 L = 2,NRNUC                                                  ! JE ADD
              ELEMNT1VEC(L) = ELEMNT1VEC(L) + DINT1VEC(IA,IB,L)*       &     ! JE ADD
                   TSHELL_R(I)                                               ! JE ADD
   21      CONTINUE                                                          ! JE ADD
           ELEMNT2 = ELEMNT2 + DINT2(IA,IB)*TSHELL_R(I)
           ELEMNT3 = ELEMNT3 + DINT3(IA,IB)*TSHELL_R(I)
           ELEMNT4 = ELEMNT4 + DINT4(IA,IB)*TSHELL_R(I)
           ELEMNT5 = ELEMNT5 + DINT5(IA,IB)*TSHELL_R(I)
           ELEMNT6 = ELEMNT6 + DINT6(IA,IB)*TSHELL_R(I)
           ELEMNT7 = ELEMNT7 + DINT7(IA,IB)*TSHELL_R(I)
        ENDDO
!
!   Return to the start of the loop
!

!*      ICI = 0
!*      DO 21 I = 1,NELMNT
!*        IRI = IROW(I)
!*         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         DO 22 J = 1,NVEC
            LOC = (J-1)*NCF
            CONTRI1 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1
            CONTRI1VEC(:) = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT1VEC(:)           ! JE ADD
            CONTRI2 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT2
            CONTRI3 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT3
            CONTRI4 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT4
            CONTRI5 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT5
            CONTRI6 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT6
            CONTRI7 = EVEC(IC+LOC)*EVEC(IR+LOC)*ELEMNT7
            IF (IR.NE.IC) THEN
               CONTRI1 = 2.0D00 * CONTRI1
               CONTRI1VEC(:) = 2.0D00 * CONTRI1VEC(:)                        ! JE ADD
               CONTRI2 = 2.0D00 * CONTRI2
               CONTRI3 = 2.0D00 * CONTRI3
               CONTRI4 = 2.0D00 * CONTRI4
               CONTRI5 = 2.0D00 * CONTRI5
               CONTRI6 = 2.0D00 * CONTRI6
               CONTRI7 = 2.0D00 * CONTRI7
            ENDIF
            DENS1(J) = DENS1(J) + CONTRI1
            DO 23 L = 2,NRNUC                                                  ! JE ADD
               DENS1VEC(J,L) = DENS1VEC(J,L) + CONTRI1VEC(L)                  ! JE ADD
   23       CONTINUE                                                          ! JE ADD
            DENS2(J) = DENS2(J) + CONTRI2
            DENS3(J) = DENS3(J) + CONTRI3
            DENS4(J) = DENS4(J) + CONTRI4
            DENS5(J) = DENS5(J) + CONTRI5
            DENS6(J) = DENS6(J) + CONTRI6
            DENS7(J) = DENS7(J) + CONTRI7
   22    CONTINUE
!*   21 CONTINUE

      GOTO 16
      ENDIF
      RETURN
      END SUBROUTINE DENSREAD_SELTZ
