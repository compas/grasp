!***********************************************************************
!                                                                      *
      SUBROUTINE SMSMCP(VINT) 
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:46:50   1/ 6/07  
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
      USE sms1_C
      USE hmat_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE gco_I 
      USE vinti_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN) :: VINT(NNNW,NNNW) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(:), pointer :: COEFF
      INTEGER, DIMENSION(:), pointer :: ICLMN, INDEX, NSWAP
      INTEGER :: NDIM, I, NWM1, KM, K, IA, NKJIA, IAP1, IB, NKJIB, KMIN, KMAX, &
         IR, IDIAG, NFILE, IOS, LAB, NCONTR, ID, IC, LOC, ICI, IRI, J 
      REAL(DOUBLE) :: PCOEFF, PICLMN, PINDEX, PNSWAP, COEF, GKAB, TEGRAL, &
         CONTRI 
      LOGICAL :: SET 
!-----------------------------------------------
!
      OPEN(UNIT=20, FILE='sms.20', FORM='UNFORMATTED', STATUS='OLD', POSITION=&
         'asis') 
 
!
!   Allocate storage that is local to this subroutine
!
      NDIM = 1 
      CALL ALLOC (COEFF, NDIM,'COEFF','SMSMCP') 
      CALL ALLOC (ICLMN, NDIM,'ICLMN','SMSMCP') 
      CALL ALLOC (INDEX, NDIM,'INDEX','SMSMCP') 
      CALL ALLOC (NSWAP, NDIM,'NSWAP','SMSMCP') 
!
!   General initializations
!
      READ (30) NELMNT 
      CALL ALLOC (IENDC, NCF + 1,'IENDC','SMSMCP') 
      CALL ALLOC (IROW, NELMNT,'IROW','SMSMCP') 
      READ (30) (IENDC(I),I=0,NCF), (IROW(I),I=1,NELMNT) 
      CLOSE(30) 
!
!   Other initializations
!
      CALL ALLOC (EMT, NELMNT,'EMT','SMSMCP') 
!
      NWM1 = NW - 1 
!
      EMT(:NELMNT) = 0.0D00 
!
!   Accumulate diagonal terms that do not require MCP coefficients
!
!                    k
!   Piece involving G (a,b) integrals
!
      KM = 0 
      K = 1 
      DO IA = 1, NWM1 
         NKJIA = NKJ(IA) 
         IAP1 = IA + 1 
         DO IB = IAP1, NW 
            NKJIB = NKJ(IB) 
            SET = .FALSE. 
            IF (NAK(IA)*NAK(IB) > 0) THEN 
               KMIN = ABS((NKJIA - NKJIB)/2) 
            ELSE 
               KMIN = ABS((NKJIA - NKJIB)/2) + 1 
            ENDIF 
            IF (MOD(K - KMIN,2) /= 0) CYCLE  
            KMAX = (NKJIA + NKJIB)/2 
            KM = MAX0(KMAX,KM) 
            IF (K<KMIN .OR. K>KMAX) CYCLE  
            DO IR = 1, NCF 
               COEF = GCO(K,IR,IA,IB) 
               IF (ABS(COEF) <= 0.0D00) CYCLE  
               IF (.NOT.SET) THEN 
!ww
                  GKAB = VINT(IA,IB)*VINT(IB,IA) 
!ww
                  SET = .TRUE. 
               ENDIF 
               IDIAG = IENDC(IR - 1) + 1 
               EMT(IDIAG) = EMT(IDIAG) - COEF*GKAB 
            END DO 
         END DO 
      END DO 
!
!   Accumulate two-electron terms that require MCP coefficients
!
      NFILE = 33 
!
      REWIND (NFILE) 
      REWIND (20) 
      READ (NFILE) 
      READ (NFILE) 
      READ (NFILE) 
!
!   The multipolarity of the integral can be deduced from the file
!   unit number
!
      K = NFILE - 32 
!
!   Attempt to read another block of data
!
   18 CONTINUE 
      READ (NFILE, IOSTAT=IOS) LAB, NCONTR 
!
      IF (IOS == 0) THEN 
!                                          k
!   Read successful; decode the labels of R (abcd)
!
         ID = MOD(LAB,KEY) 
         LAB = LAB/KEY 
         IB = MOD(LAB,KEY) 
         LAB = LAB/KEY 
         IC = MOD(LAB,KEY) 
         IA = LAB/KEY 
!
!   Compute the Vinti integrals
!
!ww
!ww            TEGRAL = VINT (IA,IB)*VINT (IC,ID)
         TEGRAL = VINT(IA,IC)*VINT(IB,ID) 
!ww
!
!   Ensure that storage is adequate to read in the rest of
!   this block
!
         IF (NCONTR > NDIM) THEN 
            CALL DALLOC (COEFF,'COEFF','SMSMCP') 
            CALL DALLOC (ICLMN,'ICLMN','SMSMCP') 
            CALL DALLOC (INDEX,'INDEX','SMSMCP') 
            CALL DALLOC (NSWAP,'NSWAP','SMSMCP') 
            NDIM = NCONTR 
            CALL ALLOC (COEFF, NDIM,'COEFF','SMSMCP') 
            CALL ALLOC (ICLMN, NDIM,'ICLMN','SMSMCP') 
            CALL ALLOC (INDEX, NDIM,'INDEX','SMSMCP') 
            CALL ALLOC (NSWAP, NDIM,'NSWAP','SMSMCP') 
         ENDIF 
!
!   Read the column index, the sparse matrix index, and the
!   coefficient for all contributions from this integral
!
         READ (NFILE) (ICLMN(I),INDEX(I),COEFF(I),I=1,NCONTR) 
         READ (20) (NSWAP(I),I=1,NCONTR) 
!
!   Store all the contributions from this integral
!
         DO I = 1, NCONTR 
            LOC = INDEX(I) 
            EMT(LOC) = EMT(LOC) - TEGRAL*COEFF(I)*(-1)**NSWAP(I) 
         END DO 
!
!   Return to the start of the loop
!
         GO TO 18 
!
      ENDIF 
!
!
!   Deallocate storage that is local to this routine
!
      CALL DALLOC (COEFF,'COEFF','SMSMCP') 
      CALL DALLOC (ICLMN,'ICLMN','SMSMCP') 
      CALL DALLOC (INDEX,'INDEX','SMSMCP') 
      CALL DALLOC (NSWAP,'NSWAP','SMSMCP') 
 
      ICI = 0 
      DO I = 1, NELMNT 
         IRI = IROW(I) 
         IF (I > IENDC(ICI)) ICI = ICI + 1 
         DO J = 1, NVEC 
            LOC = (J - 1)*NCF 
            CONTRI = EVEC(ICI + LOC)*EVEC(IRI + LOC)*EMT(I) 
            IF (IRI /= ICI) CONTRI = 2.0D00*CONTRI 
            SMSC(J) = SMSC(J) + CONTRI 
         END DO 
      END DO 
      CALL DALLOC (EMT,  'EMT',  'SMSMCP') 
      CALL DALLOC (IENDC,'IENDC','SMSMCP') 
      CALL DALLOC (IROW, 'IROW', 'SMSMCP') 
!
!   Close the angular files
!
      CLOSE(20) 
      DO I = 30, 32 + KMAXF 
         CLOSE(I) 
      END DO 
      RETURN  
      END SUBROUTINE SMSMCP 
