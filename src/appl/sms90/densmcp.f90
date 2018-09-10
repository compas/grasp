!***********************************************************************
!                                                                      *
      SUBROUTINE DENSMCP(DINT1, DINT2, DINT3, DINT4, DINT5, DINT6) 
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the  sms parameter, the electron density at the    *
!   origin and the field shift between two isotopes characterized      *
!   by fermi distributions c and a                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, ALLOC, CONVRT, DALLOC, GETYN          *
!                        ITJPO, NUCPOT, RKCO, TNSRJJ,                  *
!               [SMS92]: RINTISO, RINTDENS, VINTI                      *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:41:20   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE memory_man
      USE eigv_C
      USE mcpa_C
      USE orb_C
      USE prnt_C
      USE SMS1_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: DINT1, DINT2, &
                                                        DINT3, DINT4, &
                                                        DINT5, DINT6
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NDIM, NELMNT, I, IA, IR, IDIAG, IOS, LAB, NCONTR, IB, &
         LOC, ICI, IRI, J 
      REAL(DOUBLE), DIMENSION(:), pointer :: EMT1, EMT2, EMT3, EMT4, EMT5,  &
           EMT6, Coeff
      REAL(DOUBLE) :: QA, DIAA1, DIAA2, DIAA3, DIAA4, &
         DIAA5, DIAA6, TEGRAL1, TEGRAL2, TEGRAL3, TEGRAL4, TEGRAL5, TEGRAL6, &
         CONTRI1, CONTRI2, CONTRI3, CONTRI4, CONTRI5, CONTRI6 
      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, IENDC, IROW
      LOGICAL :: SET 
      CHARACTER(LEN=7),PARAMETER :: name='densmcp'
!-----------------------------------------------
!
!   Allocate storage that is local to this subroutine
!
      NDIM = 1 
      CALL ALLOC (COEFF, NDIM, 'COEFF', name ) 
      CALL ALLOC (ICLMN, NDIM, 'ICLMN', name) 
      CALL ALLOC (INDEX, NDIM, 'INDEX', name) 
!
      READ (30) NELMNT 
      CALL ALLOC (IENDC, 0, NCF, 'IENDC', name) 
      CALL ALLOC (IROW, NELMNT, 'IROW', name) 
      READ (30) (IENDC(I),I=0,NCF), (IROW(I),I=1,NELMNT) 
      CLOSE(30) 
!
!   Other initializations
!
      CALL ALLOC (EMT1, NELMNT, 'EMT1', name) 
      CALL ALLOC (EMT2, NELMNT, 'EMT2', name ) 
      CALL ALLOC (EMT3, NELMNT, 'EMT3', name) 
      CALL ALLOC (EMT4, NELMNT, 'EMT4', name) 
      CALL ALLOC (EMT5, NELMNT, 'EMT5', name) 
      CALL ALLOC (EMT6, NELMNT, 'EMT6', name) 
!
      EMT1(:NELMNT) = 0.0D00 
      EMT2(:NELMNT) = 0.0D00 
      EMT3(:NELMNT) = 0.0D00 
      EMT4(:NELMNT) = 0.0D00 
      EMT5(:NELMNT) = 0.0D00 
      EMT6(:NELMNT) = 0.0D00 
!
!   Accumulate diagonal terms that do not require MCP coefficients
!
!
!   Piece involving I(a,a) integrals
!
      DO IA = 1, NW 
         SET = .FALSE. 
         DO IR = 1, NCF 
            QA = DBLE(IQ(IA,IR)) 
            IF (QA <= 0.0D00) CYCLE  
            IF (.NOT.SET) THEN 
               DIAA1 = DINT1(IA,IA) 
               DIAA2 = DINT2(IA,IA) 
               DIAA3 = DINT3(IA,IA) 
               DIAA4 = DINT4(IA,IA) 
               DIAA5 = DINT5(IA,IA) 
               DIAA6 = DINT6(IA,IA) 
               SET = .TRUE. 
            ENDIF 
            IDIAG = IENDC(IR - 1) + 1 
            EMT1(IDIAG) = EMT1(IDIAG) + QA*DIAA1 
            EMT2(IDIAG) = EMT2(IDIAG) + QA*DIAA2 
            EMT3(IDIAG) = EMT3(IDIAG) + QA*DIAA3 
            EMT4(IDIAG) = EMT4(IDIAG) + QA*DIAA4 
            EMT5(IDIAG) = EMT5(IDIAG) + QA*DIAA5 
            EMT6(IDIAG) = EMT6(IDIAG) + QA*DIAA6 
         END DO 
      END DO 
!
!   Accumulate one-electron terms that require MCP coefficients
!
      REWIND (31) 
      READ (31) 
      READ (31) 
      READ (31) 
!
!   Attempt to read another block of data
!
   16 CONTINUE 
      READ (31, IOSTAT=IOS) LAB, NCONTR 
!
      IF (IOS == 0) THEN 
!
!   Read successful; decode the labels of I(ab)
!
         IA = MOD(LAB,KEY) 
         IB = LAB/KEY 
!
!   Compute I(ab)
!
         TEGRAL1 = DINT1(IA,IB) 
         TEGRAL2 = DINT2(IA,IB) 
         TEGRAL3 = DINT3(IA,IB) 
         TEGRAL4 = DINT4(IA,IB) 
         TEGRAL5 = DINT5(IA,IB) 
         TEGRAL6 = DINT6(IA,IB) 
!
!   Ensure that storage is adequate to read in the rest of
!   this block
!
         IF (NCONTR > NDIM) THEN 
            CALL DALLOC (COEFF, 'COEFF', name) 
            CALL DALLOC (ICLMN, 'ICLMN', name) 
            CALL DALLOC (INDEX, 'INDEX', name) 
            NDIM = NCONTR 
            CALL ALLOC (COEFF, NDIM, 'COEFF', name) 
            CALL ALLOC (ICLMN, NDIM, 'ICLMN', name) 
            CALL ALLOC (INDEX, NDIM, 'INDEX', name) 
         ENDIF 
!
!   Read the column index, the sparse matrix index, and the
!   coefficient for all contributions from this integral
!
         READ (31) (ICLMN(I),INDEX(I),COEFF(I),I=1,NCONTR) 
!
!   Store all the contributions from this integral
!
         DO I = 1, NCONTR 
            LOC = INDEX(I) 
            EMT1(LOC) = EMT1(LOC) + TEGRAL1*COEFF(I) 
            EMT2(LOC) = EMT2(LOC) + TEGRAL2*COEFF(I) 
            EMT3(LOC) = EMT3(LOC) + TEGRAL3*COEFF(I) 
            EMT4(LOC) = EMT4(LOC) + TEGRAL4*COEFF(I) 
            EMT5(LOC) = EMT5(LOC) + TEGRAL5*COEFF(I) 
            EMT6(LOC) = EMT6(LOC) + TEGRAL6*COEFF(I) 
         END DO 
!
!   Return to the start of the loop
!
         GO TO 16 
!
      ENDIF 
!
!   Deallocate storage that is local to this routine
!
      CALL DALLOC (COEFF, 'COEFF', name) 
      CALL DALLOC (ICLMN, 'ICLMN', name) 
      CALL DALLOC (INDEX, 'INDEX', name) 
 
      ICI = 0 
      DO I = 1, NELMNT 
         IRI = IROW(I) 
         IF (I > IENDC(ICI)) ICI = ICI + 1 
         DO J = 1, NVEC 
            LOC = (J - 1)*NCF 
            CONTRI1 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT1(I) 
            CONTRI2 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT2(I) 
            CONTRI3 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT3(I) 
            CONTRI4 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT4(I) 
            CONTRI5 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT5(I) 
            CONTRI6 = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT6(I) 
            IF (IRI /= ICI) THEN 
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
      CALL DALLOC (EMT1, 'EMT1', name) 
      CALL DALLOC (EMT2, 'EMT2', name) 
      CALL DALLOC (EMT3, 'EMT3', name) 
      CALL DALLOC (EMT4, 'EMT4', name) 
      CALL DALLOC (EMT5, 'EMT5', name) 
      CALL DALLOC (EMT6, 'EMT6', name) 
      CALL DALLOC (IENDC, 'IENDC', name) 
      CALL DALLOC (IROW, 'IROW', name) 
!
!   Close the angular files
!
      DO I = 30, 32 + KMAXF 
         CLOSE(I) 
      END DO 
 
      RETURN  
      END SUBROUTINE DENSMCP 
