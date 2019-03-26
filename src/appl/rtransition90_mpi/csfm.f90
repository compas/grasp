!***********************************************************************
!                                                                      *
      SUBROUTINE CSFM (ASFA,ASFB,LEV1,LEV2)
!                                                                      *
!   This routine calculates  the CSF Coulomb, Babuskin, and magnetic   *
!   matrix elements for  a transition  between  levels  separated by   *
!   energy OMEGA.                                                      *
!                                                                      *
!   Modified for different initial and final state orbitals            *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE biorb_C
      USE debug_C
      USE eigv_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE orb_C
      USE osc_C
      USE prnt_C
      USE syma_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE spme_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(OUT) :: asfa, asfb
      INTEGER, INTENT(IN) :: lev1, lev2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: key=KEYORB, nca=65536
      INTEGER :: ja, jb, m, nclr, ncup, ncount, nfi, ipar, jpar, &
                 ia, ib, j
      REAL(DOUBLE) :: hcoul1, hbab1, hmag1, eps, couvx, coeff, contr
      CHARACTER(LEN=4) :: jlabi, jlabj, jpari, jparj
!-----------------------------------------------------------------------
!
      EPS = 1.0D-10
!
      ASFA = 0.0D00
      ASFB = 0.0D00
      NCOUNT = 0
      NFI = 0
!
!Rasa -- start
!
!   J/pi labels for levels
!
      JLABI = LABJ(IATJPOII(LEV1))
      JLABJ = LABJ(IATJPOFF(LEV2))
      IPAR = (IASPARII(LEV1)+3)/2
      JPAR = (IASPARFF(LEV2)+3)/2
      JPARI = LABP(IPAR)
      JPARJ = LABP(JPAR)

      IF (LDBPR(18)) THEN
          WRITE (*,302) IVECFF(LEV2), JLABJ, JPARJ, IVECII(LEV1),    &
       JLABI, JPARI
          WRITE (*,301)
      ENDIF
! Rasa -- end
      DO M = 1,NINTEG
        JA = LAB(M)
        JB = MOD (JA,KEY)
        JA = JA/KEY
        NCLR = NPTR(M)+1
        NCUP = NPTR(M+1)
        IF (NCLR .GT. NCUP) GOTO 4
        IF (MOD(NKL(JA)+NKL(JB)+KK+LK,2) .NE. 0) GOTO 4
        CALL SPME (JA,JB,HCOUL1,HBAB1,HMAG1)
!
        DO J = NCLR,NCUP
!
          IA = ISLDR(J)
          IB = ISLDR1(J)
!         IB = MOD (IA,NCA)
!         IA = IA/NCA
!
!  Check on consistency
!
          IF (IA.GT.NCFII) THEN
            GOTO 131
          ELSEIF (IB.GT.NCFFF) THEN
            GOTO 131
          ENDIF
!
!  Observ that IA and IB refer to the merged list
!  whereas EVECII and EVECFF refers to the initial and
!  final state lists. Therefore we have IB -> IB-NCFII
!
          COUVX =EVECII(IA+(LEV1-1)*NCFII)*                          &
                 EVECFF(IB+(LEV2-1)*NCFFF)
!     :           EVECFF(IB-NCFII+(LEV2-1)*NCFFF)
          COEFF = XSLDR(J)
          IF (ABS(COUVX) .GT. EPS) THEN
            IF (KK .EQ. 0) THEN
              ASFA = ASFA+HCOUL1*COEFF*COUVX
              ASFB = ASFB+HBAB1*COEFF*COUVX
!Rasa -- start
              IF (LDBPR(18) .and. (dabs(HCOUL1*COEFF*COUVX)          &
                .gt. cutoff)) THEN
                  WRITE (*,300) IA, EVECII(IA+(LEV1-1)*NCFII),       &
                      IB, EVECFF(IB+(LEV2-1)*NCFFF),                 &
                      NP(JA),NH(JA),NP(JB),NH(JB),COEFF,HCOUL1,'C',  &
                      HCOUL1*COEFF*COUVX
                  WRITE(*,303) HBAB1,'B',HBAB1*COEFF*COUVX
              ENDIF
!Rasa -- end
            ELSE
              contr=HMAG1*COEFF*COUVX
              ASFA = ASFA+HMAG1*COEFF*COUVX
!Rasa -- start
              IF (LDBPR(18) .and. (abs(contr) .gt. CUTOFF)) THEN
                  WRITE (*,300) IA, EVECII(IA+(LEV1-1)*NCFII),       &
                      IB, EVECFF(IB-(LEV2-1)*NCFFF),                 &
                      NP(JA),NH(JA),NP(JB),NH(JB),COEFF,HMAG1,' ',   &
                      contr
              ENDIF
!Rasa -- end
            ENDIF
          ENDIF
 131    CONTINUE
        ENDDO
    4 ENDDO
!
      RETURN
!
  300 FORMAT (I5,2X,D12.5,I5,2X,D12.5, 2(1X,1I2,1A2),1P,D15.6,1X,D15.6, &
      1X,1A,1X,D15.6)
  301 FORMAT (4X,'i',3X,'Coeff.(i)',5X,'f',2X,'Coeff.(f)'2X,            &
      'Orb(i)',1X,'Orb(f)',2X,'MCT coeff.',4X,'Radial Factor',          &
      5X,'Contribution')
  302 FORMAT ('Upper level ',I3,A4,A4,10X,' Lower level',I3,A4,A4)
          WRITE (*,302) IVECFF(LEV2), JLABJ, JPARJ, LEV1, JLABI, JPARI
  303 FORMAT (62X,1P,D15.6,1X,1A,1X,D15.6)
  304 FORMAT(65X,'Total:')
  305 FORMAT (62X,1P,16X,1A,1X,D15.6)
!
      END
