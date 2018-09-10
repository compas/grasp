!***********************************************************************
!                                                                      *
      SUBROUTINE SPME(I, J, HCOUL, HBAB, HMAG) 
!                                                                      *
!   This routine calculates the reduced matrix elements for pair I,J   *
!   in either Coulomb/Babuskin gauge or for magnetic case.             *
!                                                                      *
!   These are defined in the  Brink and Satchler sense --- i.e. com-   *
!   patible with Pyper et al.'s MCT paper but not with Grant's radi-   *
!   ative transitions paper.                                           *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:05:40   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE bess_C
      USE biorb_C
      USE debug_C
      USE def_C
      USE grid_C
      USE orb_C
      USE osc_C
      USE tatb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE clrx_I 
      USE quad_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I,J 
      REAL(DOUBLE), INTENT(OUT) :: HCOUL, HBAB, HMAG 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NKJI, NKJJ, IPJ, NAKI, NAKJ, II, LP, LM 
      REAL(DOUBLE) :: EPS, HGAUGE, TJI, TJJ, FACT, FORM, FL, FLP, DFKI, DFKJ, &
         CLP, CLM, CIPLP, CIMLP, CIPLM, CIMLM, VALUE, CJL, CIP, TAROM, ESTHER, &
         TEKEL, PERES, TAMAR, ENOCH, SETH, SHEM, DARIUS, CYRUS, DANIEL, BOAZ, &
         ESAU, AARON, GIMEL 
!-----------------------------------------------
 
!
      EPS = 1.0D-10 
!
      IF (LDBPR(12) .OR. LDBPR(13) .OR. LDBPR(14) .OR. LDBPR(15)) WRITE (99, &
         303) I, J 
!
      HMAG = 0.0D00 
      HCOUL = 0.0D00 
      HGAUGE = 0.0D00 
      HBAB = 0.0D00 
!
!  Evaluate factor multiplying Mbar(a,b)
!
      NKJI = NKJ(I) 
      NKJJ = NKJ(J) 
      TJI = DBLE(NKJI) 
      TJJ = DBLE(NKJJ) 
      FACT = 1.0D00 
      IPJ = (NKJI + 1)/2 + NKJJ 
      IF (MOD(IPJ,2) /= 0) FACT = -FACT 
      NAKI = NAK(I) 
      NAKJ = NAK(J) 
      FORM = FACT*SQRT(TJJ + 1.0D00)*CLRX(NAKJ,LK,NAKI) 
      FL = DBLE(LK) 
      FLP = FL + 1.0D00 
      DFKI = DBLE(NAKI) 
      DFKJ = DBLE(NAKJ) 
!
      MTP = MIN(MFII(NNII(I)),MFFF(NNFF(J))) 
!
      IF (KK == 0) THEN 
!
         IF (ABS(FORM) > EPS) THEN 
!
!  To pick the right initial and final state radial functions convert
!  from orbital order of the merged list to the orbital orders of the
!  initial and final state lists, respectively.
!  N(I) -> NNII(I), N(J) -> NNFF(J)
!
            DO II = 1, N 
               TB(II) = PFII(II,NNII(I))*QFFF(II,NNFF(J)) + QFII(II,NNII(I))*&
                  PFFF(II,NNFF(J)) 
               TC(II) = PFII(II,NNII(I))*QFFF(II,NNFF(J)) - QFII(II,NNII(I))*&
                  PFFF(II,NNFF(J)) 
! zou
               TB(II) = TB(II)*C/CVAC 
               TC(II) = TC(II)*C/CVAC 
               TD(II) = PFII(II,NNII(I))*PFFF(II,NNFF(J)) + QFII(II,NNII(I))*&
                  QFFF(II,NNFF(J))*(C/CVAC)**2 
!           TD(II) = PFII(II,NNII(I))*PFFF(II,NNFF(J)) +
!    :               QFII(II,NNII(I))*QFFF(II,NNFF(J))
! zou
            END DO 
            LP = LK + 1 
            LM = LK - 1 
!
!  Calculate Coulomb coefficients
!
            CLP = SQRT(FL/FLP) 
            CLM = -SQRT(FLP/FL) 
            CIPLP = CLP*(DFKI - DFKJ) 
            CIMLP = CLP*FLP 
            CIPLM = CLM*(DFKI - DFKJ) 
            CIMLM = -CLM*FL 
!
!  Tabulate Coulomb integrand
!
            TA(1) = 0.0D00 
            TA(2:MTP) = RP(2:MTP)*(BJ(2:MTP,3)*(CIPLP*TB(2:MTP)+CIMLP*TC(2:MTP)&
               )+BJ(2:MTP,1)*(CIPLM*TB(2:MTP)+CIMLM*TC(2:MTP))) 
            CALL QUAD (VALUE) 
            HCOUL = FORM*VALUE 
!
!  Calculate gauge dependent coefficients
!
            CJL = -(FL + FLP) 
            CIP = DFKI - DFKJ 
!
!  Tabulate gauge dependent integrand
!
            TA(1) = 0.0D00 
            TA(2:MTP) = RP(2:MTP)*(BJ(2:MTP,2)*(CJL*TD(2:MTP))+BJ(2:MTP,3)*(CIP&
               *TB(2:MTP)+FLP*TC(2:MTP))+BJ(2:MTP,1)*(CIP*TB(2:MTP)-FL*TC(2:MTP&
               ))) 
!
!  Print gauge dependent integrand if requested
!
            IF (LDBPR(13)) WRITE (99, 300) I, J, (II,TA(II),II=1,N) 
            CALL QUAD (VALUE) 
            HGAUGE = FORM*VALUE 
            HBAB = HCOUL + SQRT(FLP/FL)*HGAUGE 
!
!  Print Coulomb and gauge dependent integrals if requested
!
            IF (LDBPR(12)) WRITE (99, 301) I, J, HCOUL, HGAUGE, HBAB 
!
!  Calculate the contibutions from various terms if requested
!
            IF (LDBPR(14) .OR. LDBPR(15)) THEN 
!
               TA(1) = 0.0D00 
               TA(2:MTP) = TD(2:MTP)*BJ(2:MTP,2)*RP(2:MTP) 
               IF (LDBPR(15)) WRITE (99, 304) (II,TA(II),II=1,N) 
               CALL QUAD (TAROM) 
               IF (LDBPR(14)) WRITE (99, 305) TAROM 
               ESTHER = -CJL 
               IF (LDBPR(14)) WRITE (99, 306) ESTHER 
               TEKEL = ESTHER*TAROM 
               IF (LDBPR(14)) WRITE (99, 307) TEKEL 
!
               TA(1) = 0.0D00 
               TA(2:MTP) = -TB(2:MTP)*BJ(2:MTP,3)*RP(2:MTP) 
               IF (LDBPR(15)) WRITE (99, 304) (II,TA(II),II=1,N) 
               CALL QUAD (PERES) 
               IF (LDBPR(14)) THEN 
                  WRITE (99, 308) PERES 
                  WRITE (99, 306) CIP 
               ENDIF 
               TAMAR = CIP*PERES 
               IF (LDBPR(14)) WRITE (99, 307) TAMAR 
!
               TA(1) = 0.0D00 
               TA(2:MTP) = -TB(2:MTP)*BJ(2:MTP,1)*RP(2:MTP) 
               IF (LDBPR(15)) WRITE (99, 304) (II,TA(II),II=1,N) 
               CALL QUAD (ENOCH) 
               IF (LDBPR(14)) THEN 
                  WRITE (99, 309) ENOCH 
                  WRITE (99, 306) CIP 
               ENDIF 
               SETH = CIP*ENOCH 
               IF (LDBPR(14)) WRITE (99, 307) SETH 
!
               TA(1) = 0.0D00 
               TA(2:MTP) = BJ(2:MTP,1)*TC(2:MTP)*RP(2:MTP) 
               IF (LDBPR(15)) WRITE (99, 304) (II,TA(II),II=1,N) 
               CALL QUAD (SHEM) 
               IF (LDBPR(14)) THEN 
                  WRITE (99, 310) SHEM 
                  WRITE (99, 306) FL 
               ENDIF 
               DARIUS = FL*SHEM 
               IF (LDBPR(14)) WRITE (99, 307) DARIUS 
!
               TA(1) = 0.0D00 
               TA(2:MTP) = -BJ(2:MTP,3)*TC(2:MTP)*RP(2:MTP) 
               IF (LDBPR(15)) WRITE (99, 304) (II,TA(II),II=1,N) 
               CALL QUAD (CYRUS) 
               IF (LDBPR(14)) WRITE (99, 311) CYRUS 
               DANIEL = FLP 
               IF (LDBPR(14)) WRITE (99, 306) DANIEL 
               BOAZ = DANIEL*CYRUS 
               IF (LDBPR(14)) WRITE (99, 307) BOAZ 
!
               ESAU = TAMAR + SETH + DARIUS + BOAZ 
               IF (LDBPR(14)) WRITE (99, 312) ESAU, TEKEL 
               AARON = ESAU + TEKEL 
               GIMEL = -VALUE 
               IF (LDBPR(14)) THEN 
                  WRITE (99, 313) AARON, GIMEL 
                  WRITE (99, 314) 
                  WRITE (99, 315) PERES, CIP, TAMAR 
                  WRITE (99, 316) ENOCH, CIP, SETH 
                  WRITE (99, 317) SHEM, FL, DARIUS 
                  WRITE (99, 318) CYRUS, DANIEL, BOAZ 
                  WRITE (99, 319) ESAU 
                  WRITE (99, 320) TAROM, ESTHER, TEKEL 
                  WRITE (99, 321) GIMEL, AARON 
               ENDIF 
!
            ENDIF 
!
         ENDIF 
!
         IF (LDBPR(12)) WRITE (99, 301) I, J, HCOUL, HGAUGE, HBAB 
!
      ELSE 
!
!   Tabulate magnetic integrand
!
         IF (ABS(FORM) > EPS) THEN 
!
            TA(1) = 0.0D00 
            DO II = 2, MTP 
!zou
               TA(II) = (PFII(II,NNII(I))*QFFF(II,NNFF(J))+QFII(II,NNII(I))*&
                  PFFF(II,NNFF(J)))*BJ(II,2)*RP(II)*C/CVAC 
!           TA(II) =  (PFII(II,NNII(I))*QFFF(II,NNFF(J)) +
!    :      QFII(II,NNII(I))*PFFF(II,NNFF(J)))*BJ(II,2)*RP(II)
!zou
            END DO 
            CALL QUAD (VALUE) 
            IF (LDBPR(14)) WRITE (99, 322) VALUE 
            HMAG = -VALUE*(FL + FLP)*(DFKI + DFKJ)*FORM/SQRT(FL*FLP) 
!
         ENDIF 
!
         IF (LDBPR(12)) WRITE (99, 302) HMAG 
!
      ENDIF 
!
      RETURN  
!
  300 FORMAT(/,1X,20X,'Local form of gauge dependent integral for',' orbitals',&
         I3,' and',I3,/,/,100(1X,7(I3,2X,1P,D11.3,2X),/)) 
  301 FORMAT(/,1X,'Orbital pair (',I2,',',I2,') integrals:','Coul. gauge ',1P,D&
         13.3,3X,'Gauge contribution',D13.3,3X,'Bab.gauge',D13.3) 
  302 FORMAT(/,1X,'Magnetic single particle matrix element ',1P,D13.3) 
  303 FORMAT(/,1X,30('+'),'   SPME called for orbitals ',I4,'   and',I4,3X,30(&
         '+')) 
  304 FORMAT(/,100(1X,7(I3,2X,1P,D11.3,2X),/),/) 
  305 FORMAT(/,' Integral J(L)  =  ',1P,D28.12) 
  306 FORMAT(' Factor multiplying this ',1P,D20.12) 
  307 FORMAT(' Resultant contribution  ',1P,D20.12,/) 
  308 FORMAT(' I + (L+1) integral  = ',1P,D24.12) 
  309 FORMAT(' I + (L-1) integral  = ',1P,D24.12) 
  310 FORMAT(' I - (L-1) integral  = ',1P,D24.12) 
  311 FORMAT(' I - (L+1) integral  = ',1P,D24.12) 
  312 FORMAT(' These last four add up to ',1P,D20.12,&
         '    compared to the first which is ',D20.12) 
  313 FORMAT(' These all add up to ',1P,D20.12,' which should equal ',D20.12) 
  314 FORMAT(/,18X,'Integral',11X,'Coefficient',8X,'Contribution') 
  315 FORMAT(' I + (L+1)',1P,3D20.10) 
  316 FORMAT(' I + (L-1)',1P,3D20.10) 
  317 FORMAT(' I - (L-1)',1P,3D20.10) 
  318 FORMAT(' I - (L+1)',1P,3D20.10) 
  319 FORMAT(54X,16('-'),/,50X,1P,D20.10) 
  320 FORMAT(' J (L)',4X,1P,3D20.10) 
  321 FORMAT(54X,16('-'),/,' Correct figure',1P,D20.10,15X,D20.10) 
  322 FORMAT(' I + (L)  integral',1P,D20.10) 
      RETURN  
!
      END SUBROUTINE SPME 
