!***********************************************************************
!                                                                      *
      SUBROUTINE HFSGG
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the  hyperfine structure parameters.               *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, DRACAH, ISPAR, ITJPO,  *
!                        TNSRJJ.                                       *
!               [HFS92]: MATELT, RINT, RINTHF.                         *
!                                                                      *
!   Written by Per Jonsson and Farid A. Parpia                         *
!                                                                      *
!   Modified by Per Jonsson to evaluate g_j factors                    *
!                                                                      *
!                                         Last revision: 22 Oct 1999   *
!   Modified by Gediminas Gaigalas for new spin-angular integration.   *
!                                         Last revision:    Nov 2017   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:35:13   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW
      USE memory_man
      USE decide_C
      USE DEF_C,           ONLY: AUMAMU, CVAC, EMPAM, RBCM, AUCM,     &
                                 CCMS, B1=>AUMAMU
      USE EIGV_C 
      USE foparm_C
      USE jlabl_C,               LABJ=>JLBR, LABP=>JLBP
      USE nsmdat_C,        ONLY: SQN, DMOMNM, QMOMB,                  &
                                 HFSI=>SQN, HFSD=>DMOMNM, HFSQ=>QMOMB
      USE orb_C
      USE prnt_C
      USE OPT6_C 
      USE syma_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rinthf_I 
      USE rint_I 
      USE matelt_I 
      USE convrt_I 
      USE ispar_I 
      USE itjpo_I 
      USE oneparticlejj_I 
      USE gracah1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-10 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: FFMIN, FFMAX, FF, I, J, KT, IPT, IC, LCNUM, IR, ISPARC, ITJPOC&
         , ITJPOR, IDIFF, IA, IB, K, KK, LOC1, LOC2, II, JJ, JJII, JB, JA, JJB&
         , JJA, IFLAG 
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL 
      REAL(DOUBLE), DIMENSION(2,NNNW,NNNW) :: RINTME, AMELT 
      REAL(DOUBLE), DIMENSION(NNNW,NNNW) :: RINTGJ, RINTDGJ, GJMELT, DGJMELT 
!     .. Local pointer arrays
      REAL(DOUBLE), DIMENSION(:,:), pointer :: HFC
      REAL(DOUBLE),  DIMENSION(:), pointer :: GJC, DGJC
      REAL(DOUBLE) ::  APART, GJPART, DGJPART, ELEMNT,  ELEMNTGJ, ELEMNTDGJ,&
         CONTR, CONTRGJ, CONTRDGJ, AUMHZ, BARNAU, DNMAU, GFAC, HFAC, FJ, &
         GJA1, AFA1, AFA2, BFA1, BFA2, BFA3, GJ, DGJ, TILDE1, &
         TILDE2, FACTOR1, FACTOR2, RAC1, RAC2, HFSELT1, HFSELT2 
      CHARACTER :: CNUM*11 
!-----------------------------------------------
!
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!
!   Allocate storage for local arrays
!
      CALL ALLOC (HFC, 5, NVEC*NVEC, 'HFC', 'HFS') 
      CALL ALLOC (GJC, NVEC*NVEC, 'GJC', 'HFS') 
      CALL ALLOC (DGJC, NVEC*NVEC, 'DGLC', 'HFS') 
!
!   Initialise
!
      HFC(:,:NVEC*NVEC) = 0.0D00 
!
      GJC(:NVEC*NVEC) = 0.0D00 
      DGJC(:NVEC*NVEC) = 0.0D00 
 
!
!   Calculate and save the radial integrals and angular
!   matrix elements for the two multipolarities
!
      DO KT = 1, 2 
         DO I = 1, NW 
            DO J = 1, NW 
               IF (KT == 1) THEN 
                  RINTME(KT,I,J) = RINTHF(I,J,-2) 
                  RINTGJ(I,J) = RINTHF(I,J,1) 
               ELSE 
                  RINTME(KT,I,J) = RINT(I,J,-3) 
                  RINTDGJ(I,J) = RINT(I,J,0) 
               ENDIF 
               CALL MATELT (I, KT, J, APART, GJPART, DGJPART) 
               AMELT(KT,I,J) = APART 
               IF (KT /= 1) CYCLE  
               GJMELT(I,J) = GJPART 
               DGJMELT(I,J) = DGJPART 
            END DO 
         END DO 
      END DO 
!
!   Set the parity of the one-body operators
!
      IPT = 1 
!
!   Sweep through the Hamiltonian matrix to determine the
!   diagonal and off-diagonal hyperfine constants
!
      DO IC = 1, NCF 
!
!   Output IC on the screen to show how far the calculation has preceede
!
         IF (MOD(IC,100) == 0) THEN 
            CALL CONVRT (IC, CNUM, LCNUM) 
            WRITE (6, *) 'Column '//CNUM(1:LCNUM)//' complete;' 
         ENDIF 
!
         DO IR = 1, NCF 
!
!   If LFORDR is .TRUE., a `first order' calculation is indicated;
!   only the CSFs with serial numbers exceeding IC are treated specially
!   only diagonal elements are evaluated for the `first order' CSFs
!
            IF (LFORDR .AND. IC>ICCUT .AND. IC/=IR) CYCLE  
!
            ISPARC = ISPAR(IC) 
            ITJPOC = ITJPO(IC) 
            ITJPOR = ITJPO(IR) 
            IDIFF = ITJPOC - ITJPOR 
!
!   Loop over the multipolarities
!
            DO KT = 1, 2 
!
!   Initialise the accumulator
!
               ELEMNT = 0.0D00 
               ELEMNTGJ = 0.0D00 
               ELEMNTDGJ = 0.0D00 
!
!   Consider 3 cases
!               (k)
!   (1) < J || T   || J >    , k = 1,2  and IR >= IC
!
!               (k)
!   (2) < J || T   || J-1 >  , k = 1,2
!
!               (k)
!   (3) < J || T   || J-2 >  , k = 2
!
               IF (.NOT.(IDIFF==0 .AND. IR>=IC .OR. IDIFF==2 .OR. IDIFF==4&
                   .AND. KT==2)) CYCLE  
!
                CALL ONEPARTICLEJJ(KT,IPT,IC,IR,IA,IB,TSHELL)
!GG               CALL TNSRJJ (KT, IPT, IC, IR, IA, IB, TSHELL) 
!
!   Accumulate the contribution from the one-body operators;
!
               IF (IA /= 0) THEN 
                  IF (IA == IB) THEN 
                     IF (KT/=1 .OR. IDIFF/=0) THEN 
                        DO IA = 1, NW 
                           IF (ABS(TSHELL(IA)) <= CUTOFF) CYCLE  
                           ELEMNT = ELEMNT + AMELT(KT,IA,IA)*RINTME(KT,IA,IA)*&
                              TSHELL(IA) 
                           CYCLE  
                        END DO 
                     ELSE 
                        DO IA = 1, NW 
                           IF (ABS(TSHELL(IA)) <= CUTOFF) CYCLE  
                           ELEMNT = ELEMNT + AMELT(KT,IA,IA)*RINTME(KT,IA,IA)*&
                              TSHELL(IA) 
                           ELEMNTGJ = ELEMNTGJ + GJMELT(IA,IA)*RINTGJ(IA,IA)*&
                              TSHELL(IA) 
                           ELEMNTDGJ = ELEMNTDGJ + DGJMELT(IA,IA)*RINTDGJ(IA,IA&
                              )*TSHELL(IA) 
                        END DO 
                     ENDIF 
                  ELSE 
                     IF (ABS(TSHELL(1)) > CUTOFF) THEN 
                        ELEMNT = ELEMNT + AMELT(KT,IA,IB)*RINTME(KT,IA,IB)*&
                           TSHELL(1) 
                        IF (KT==1 .AND. IDIFF==0) THEN 
                           ELEMNTGJ = ELEMNTGJ + GJMELT(IA,IB)*RINTGJ(IA,IB)*&
                              TSHELL(1) 
                           ELEMNTDGJ = ELEMNTDGJ + DGJMELT(IA,IB)*RINTDGJ(IA,IB&
                              )*TSHELL(1) 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
!
!   Multiply with the configuration expansion coefficients and add the
!   contributions from the matrix elements to obtain total contributions
!
               DO K = 1, NVEC 
                  DO KK = 1, NVEC 
                     LOC1 = (K - 1)*NCF 
                     LOC2 = (KK - 1)*NCF 
                     IF (IDIFF==0 .AND. IR/=IC) THEN 
                        CONTR = ELEMNT*(EVEC(IC + LOC1)*EVEC(IR + LOC2) + EVEC(&
                           IR + LOC1)*EVEC(IC + LOC2)) 
                        CONTRGJ = ELEMNTGJ*(EVEC(IC + LOC1)*EVEC(IR + LOC2) + &
                           EVEC(IR + LOC1)*EVEC(IC + LOC2)) 
                        CONTRDGJ = ELEMNTDGJ*(EVEC(IC + LOC1)*EVEC(IR + LOC2)&
                            + EVEC(IR + LOC1)*EVEC(IC + LOC2)) 
                     ELSE 
                        CONTR = ELEMNT*EVEC(IC + LOC1)*EVEC(IR + LOC2) 
                        CONTRGJ = ELEMNTGJ*EVEC(IC + LOC1)*EVEC(IR + LOC2) 
                        CONTRDGJ = ELEMNTDGJ*EVEC(IC + LOC1)*EVEC(IR + LOC2) 
                     ENDIF 
!
!   Magnetic dipole and the two operators of the g_j factor
!
                     IF (KT == 1) THEN 
                        IF (IDIFF == 0) THEN 
                           HFC(1,NVEC*(K-1)+KK) = HFC(1,NVEC*(K-1)+KK) + CONTR 
                           GJC(NVEC*(K-1)+KK) = GJC(NVEC*(K-1)+KK) + CONTRGJ 
                           DGJC(NVEC*(K-1)+KK) = DGJC(NVEC*(K-1)+KK) + CONTRDGJ 
                        ELSE IF (ITJPOC - ITJPOR == 2) THEN 
                           HFC(2,NVEC*(K-1)+KK) = HFC(2,NVEC*(K-1)+KK) + CONTR 
                        ENDIF 
!
!   Electric quadrupole
!
                     ELSE IF (KT == 2) THEN 
                        SELECT CASE (IDIFF)  
                        CASE (0)  
                           HFC(3,NVEC*(K-1)+KK) = HFC(3,NVEC*(K-1)+KK) + CONTR 
                        CASE (2)  
                           HFC(4,NVEC*(K-1)+KK) = HFC(4,NVEC*(K-1)+KK) + CONTR 
                        CASE (4)  
                           HFC(5,NVEC*(K-1)+KK) = HFC(5,NVEC*(K-1)+KK) + CONTR 
                        END SELECT 
                     ENDIF 
                  END DO 
               END DO 
!
            END DO 
!
         END DO 
      END DO 
!
!   These are the conversion factors to obtain the hyperfine
!   constants in MHz
!
      AUMHZ = AUCM*CCMS*1.0D-06 
      BARNAU = 1.0D-24/RBCM**2 
      DNMAU = B1/(2.0D00*CVAC*EMPAM) 
!
      GFAC = AUMHZ*DNMAU*HFSD/HFSI 
      HFAC = AUMHZ*2.0D00*HFSQ*BARNAU 
!
!   Output the hyperfine interaction constants
!
      WRITE (24, 302) 
      WRITE (29, 402) 
!
      DO I = 1, NVEC 
         DO II = 1, NVEC 
!
            JJ = IATJPO(I) 
            JJII = IATJPO(II) 
!
            IF (.NOT.(JJ==JJII .AND. JJII>1 .OR. JJ>JJII)) CYCLE  
!
            FJ = 0.5D00*DBLE(JJ - 1) 
!
            GJA1 = SQRT(1.0D00/(FJ*(FJ + 1.0D00))) 
 
            AFA1 = GFAC*GJA1 
 
            IF (JJ == 2) THEN 
               AFA2 = 0.D000 
            ELSE 
               AFA2 = GFAC*SQRT(1.0D00/(FJ*(2.0D00*FJ - 1.0D00))) 
            ENDIF 
            BFA1 = HFAC*SQRT((FJ*(2.0D00*FJ - 1.0D00))/((FJ + 1.0D00)*(2.0D00*&
               FJ + 3.0D00))) 
            IF (JJ == 2) THEN 
               BFA2 = 0.0D00 
            ELSE 
               BFA2 = 0.25D00*HFAC*SQRT((FJ*(FJ - 1.0D00))/((FJ + 1.0D00)*(&
                  2.0D00*FJ - 1.0D00))) 
            ENDIF 
            IF (JJ == 4) THEN 
               BFA3 = 0.0D00 
            ELSE 
               BFA3 = 0.125D00*HFAC*SQRT((FJ*(FJ - 1.0D000)*(2.0D00*FJ - 1.0D00&
                  ))/(2.0D00*FJ - 3.0D00)) 
            ENDIF 
!
!   Diagonal (J,J) A and B factors
!
            SELECT CASE (JJ - JJII)  
            CASE (0)  
               IF (I <= II) THEN 
!                GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
!                DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)

                  WRITE (24, 303) IVEC(I), LABJ(JJ), LABP((IASPAR(I)+3)/2), &
                     IVEC(II), LABJ(JJII), LABP((IASPAR(II)+3)/2), AFA1*HFC(1,&
                     NVEC*(I-1)+II), BFA1*HFC(3,NVEC*(I-1)+II) 
!
!   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
!
                  IF (I == II) THEN 
                     GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II) 
                     DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II) 
                     WRITE (29, 403) IVEC(I), LABJ(JJ), LABP((IASPAR(I)+3)/2), &
                        AFA1*HFC(1,NVEC*(I-1)+II), BFA1*HFC(3,NVEC*(I-1)+II), &
                        GJ, DGJ, GJ + DGJ 
                  ENDIF 
               ENDIF 
!
!   Off diagonal (J,J-1) A and B factors
!
            CASE (2)  
               WRITE (24, 303) IVEC(I), LABJ(JJ), LABP((IASPAR(I)+3)/2), IVEC(&
                  II), LABJ(JJII), LABP((IASPAR(II)+3)/2), AFA2*HFC(2,NVEC*(I-1&
                  )+II), BFA2*HFC(4,NVEC*(I-1)+II) 
!
!   Off diagonal (J,J-2) B factor
!
            CASE (4)  
               WRITE (24, 303) IVEC(I), LABJ(JJ), LABP((IASPAR(I)+3)/2), IVEC(&
                  II), LABJ(JJII), LABP((IASPAR(II)+3)/2), 0.0D00, BFA3*HFC(5,&
                  NVEC*(I-1)+II) 
            END SELECT 
!
         END DO 
      END DO 
!
!   These are the factors needed to obtain the F-dependent hyperfine
!   matrix elements in Hartrees
!
      TILDE1 = SQRT((HFSI + 1.0D00)/HFSI)*HFSD*0.5D00*B1/(EMPAM*CVAC) 
      IF (HFSI > 0.6D00) THEN 
         TILDE2 = SQRT((HFSI + 1.0D00)*(2.0D00*HFSI + 3.D0)/((2.0D00*HFSI - &
            1.0D00)*HFSI)) 
         TILDE2 = TILDE2*HFSQ*0.5D00*BARNAU 
      ELSE 
         TILDE2 = 0.0D00 
      ENDIF 
!
!     II = 2*I
!
      II = NINT(2.0D00*HFSI) 
!
!   Calculate the F-dependent matrix elemnts
!
!   Loop over the states
!
      DO JB = 1, NVEC 
         DO JA = 1, NVEC 
            JJB = IATJPO(JB) - 1 
            JJA = IATJPO(JA) - 1 
            IF (.NOT.(JJA==JJB .AND. JJA>0 .OR. JJB>JJA)) CYCLE  
!
!   Determine the possible F quantum numbers for the matrix element
!
            FFMIN = MAX(ABS(JJA - II),ABS(JJB - II)) 
            FFMAX = MIN(JJA + II,JJB + II) 
!
!   Loop over the possible F quantum numbers
!
            IFLAG = 0 
            DO FF = FFMIN, FFMAX, 2 
!
!   Phase factor
!
               IF (MOD((II + JJA - FF)/2,2) == 1) THEN 
                  FACTOR1 = -1.0D00 
               ELSE 
                  FACTOR1 = 1.0D00 
               ENDIF 
               FACTOR2 = FACTOR1*SQRT((DBLE(JJB) + 1.0D00)*(DBLE(II) + 1.0D00)) 
!
!   Determine the Racah W coefficients.
!
!GG               CALL DRACAH (II, JJA, II, JJB, FF, 2, RAC1) 
              CALL GRACAH1 (II,JJA,II,JJB,FF,2,RAC1)
!GG               CALL DRACAH (II, JJA, II, JJB, FF, 4, RAC2) 
              CALL GRACAH1 (II,JJA,II,JJB,FF,4,RAC2)
!
!   Obtain and output matrix elements for J,J
!
               IF (JJA - JJB == 0) THEN 
                  HFSELT1 = FACTOR2*RAC1*HFC(1,NVEC*(JB-1)+JA)*TILDE1 
                  HFSELT2 = FACTOR2*RAC2*HFC(3,NVEC*(JB-1)+JA)*TILDE2 
                  IF (ABS(HFSELT1 + HFSELT2) > CUTOFF*10.0D-05) THEN 
                     IF (IFLAG == 0) WRITE (24, 304) 
                     IFLAG = 1 
                     WRITE (24, 305) IVEC(JB), LABJ(JJB+1), LABP((IASPAR(JB)+3)&
                        /2), IVEC(JA), LABJ(JJA+1), LABP((IASPAR(JA)+3)/2), &
                        LABJ(FF+1), HFSELT1 + HFSELT2 
                  ENDIF 
!
!   Obtain and output matrix elements for J,J-1
!
               ELSE IF (ABS(JJA - JJB) == 2) THEN 
                  HFSELT1 = FACTOR2*RAC1*HFC(2,NVEC*(JB-1)+JA)*TILDE1 
                  HFSELT2 = FACTOR2*RAC2*HFC(4,NVEC*(JB-1)+JA)*TILDE2 
                  IF (ABS(HFSELT1 + HFSELT2) > CUTOFF*10.0D-05) THEN 
                     IF (IFLAG == 0) WRITE (24, 304) 
                     IFLAG = 1 
                     WRITE (24, 305) IVEC(JB), LABJ(JJB+1), LABP((IASPAR(JB)+3)&
                        /2), IVEC(JA), LABJ(JJA+1), LABP((IASPAR(JA)+3)/2), &
                        LABJ(FF+1), HFSELT1 + HFSELT2 
                  ENDIF 
!
!   Obtain and output matrix elements for J,J-2
!
               ELSE IF (ABS(JJA - JJB) == 4) THEN 
                  HFSELT1 = 0.0D00 
                  HFSELT2 = FACTOR2*RAC2*HFC(5,NVEC*(JB-1)+JA)*TILDE2 
                  IF (ABS(HFSELT1 + HFSELT2) > CUTOFF*10.0D-05) THEN 
                     IF (IFLAG == 0) WRITE (24, 304) 
                     IFLAG = 1 
                     WRITE (24, 305) IVEC(JB), LABJ(JJB+1), LABP((IASPAR(JB)+3)&
                        /2), IVEC(JA), LABJ(JJA+1), LABP((IASPAR(JA)+3)/2), &
                        LABJ(FF+1), HFSELT1 + HFSELT2 
                  ENDIF 
               ELSE 
                  HFSELT1 = 0.0D00 
                  HFSELT2 = 0.0D00 
               ENDIF 
            END DO 
         END DO 
      END DO 
!
      CALL DALLOC (HFC, 'HFC', 'HFS') 
      RETURN  
!
  302 FORMAT(/,/,' Interaction constants:'/,/,&
         ' Level1  J Parity  Level2  J Parity',8X,'A (MHz)',13X,'B (MHz)'/) 
  402 FORMAT(/,/,' Interaction constants:'/,/,' Level1  J Parity ',8X,'A (MHz)'&
         ,13X,'B (MHz)',13X,'g_J',14X,'delta g_J',11X,'total g_J'/) 
  303 FORMAT(1X,1I3,5X,2A4,2X,1I3,5X,2A4,1P,2D20.10) 
  403 FORMAT(1X,1I3,5X,2A4,1P,5D20.10) 
  304 FORMAT(/,/,' Matrix elements:'/,/,&
         ' Level1  J Parity  Level2  J Parity    F ',4X,'Matrix element (a.u.)'&
         /) 
  305 FORMAT(1X,1I3,5X,2A4,2X,1I3,5X,2A4,1X,A4,4X,1P,1D20.10) 
      RETURN  
!
      END SUBROUTINE HFSGG
