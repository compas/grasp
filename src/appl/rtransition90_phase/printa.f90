!***********************************************************************
!                                                                      *
      SUBROUTINE PRINTA(ASFA, ASFB, I, J, OMEGA, FACTOR, LINES, LSAME)
!                                                                      *
!   This  routine  prints the basic oscillator strength  information   *
!   for transitions between level I and level J.                       *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:39:37   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE def_C
      USE CUTO_C
      USE JLABL_C, LABJ=>JLBR, LABP=> JLBP
      USE OSC_C
      USE SYMA_C
      USE PRNT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN) :: LSAME
      INTEGER  :: I, J
      INTEGER, INTENT(INOUT) :: LINES
      REAL(DOUBLE), INTENT(IN) :: ASFA, ASFB
      REAL(DOUBLE), INTENT(IN) :: OMEGA, FACTOR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IPAR, JPAR, ISIGNA, ISIGNB
      REAL(DOUBLE), DIMENSION(10) :: DBLFAC
      REAL(DOUBLE) :: STFAC, DLL1, OMC, FAAU, FOSC, FBAU, ENG, GGFACTOR, &
         ACSQ, ABSQ, AC, AB, BC, BB, OSCC, OSCB, SA, SB, AMS, AM, BM, OSCM
      CHARACTER(LEN=4) :: JLABI, JLABJ, JPARI, JPARJ
      CHARACTER(LEN=2) :: F1,F2
!-----------------------------------------------
!
!
!  *****  DBLFAC(I) = (2I+1)!!
!
      DATA DBLFAC/ 3.0000000000D00, 1.5000000000D01, 1.0500000000D02, &
         9.4500000000D02, 1.0395000000D04, 1.3513500000D05, 2.0270250000D06, &
         3.4459425000D07, 6.5472907500D08, 1.3749310575D10/
!
!   Evaluate statistical factors and constants
!
      STFAC = IATJPOII(I)
      DLL1 = DBLE(LK + LK + 1)
      STFAC = STFAC/DLL1
!
      OMC = OMEGA/CVAC
!     OMC = OMEGA/C
      FAAU = 2.0D00*OMC*STFAC/DBLE(IATJPOFF(J))
      FOSC = CVAC*STFAC/OMC
!     FOSC = C*STFAC/OMC
      FBAU = PI*FOSC/OMEGA
      ENG = OMEGA*FACTOR
      IF (LTC(1)) ENG = FACTOR/OMEGA
!
!   J/pi labels for levels
!
      JLABI = LABJ(IATJPOII(I))
      JLABJ = LABJ(IATJPOFF(J))
      IPAR = (IASPARII(I) + 3)/2
      JPAR = (IASPARFF(J) + 3)/2
      JPARI = LABP(IPAR)
      JPARJ = LABP(JPAR)
!GG
      GGFACTOR = CVAC**(2*LK - 2)*DBLE(LK)*DBLFAC(LK)**2/((2.0D00*DBLE(LK) + &
         1.0D00)*(DBLE(LK) + 1.0D00)*ABS(OMEGA)**(2*LK - 1))
!GG-end
!
!   Calculate Einstein A and B coefficients and oscillator strengths
!
      IF (KK == 0) THEN
!
!   Electric multipoles
!
!   In atomic units
!
         IF (ASFA < 0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         IF (ASFB < 0) THEN
            ISIGNB = -1
         ELSE
            ISIGNB = 1
         ENDIF

         ACSQ = ASFA**2
         ABSQ = ASFB**2
         AC = ACSQ*FAAU
         AB = ABSQ*FAAU
         BC = ACSQ*FBAU
         BB = ABSQ*FBAU
         OSCC = ACSQ*FOSC
         OSCB = ABSQ*FOSC
         SA = OSCC*GGFACTOR
         SB = OSCB*GGFACTOR
!
!   Convert to SI units if option 5 not set
!
         IF (.NOT.LTC(7)) THEN
            AC = AC*FASI
            AB = AB*FASI
            BC = BC*FBSI
            BB = BB*FBSI
         ENDIF
!
!   Accumulate total of A coefficients
!
         TOTC(J) = TOTC(J) + AC
         TOTB(J) = TOTB(J) + AB
!
!   Print information if both AC and AB are greater than CUTOFF
!
         IF (ABS(AC)>=CUTOFF .AND. ABS(AB)>=CUTOFF) THEN
            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f2'
                  F2 = 'f1'
               END IF
!
            WRITE (24, 300) F1,IVECFF(J),JLABJ,JPARJ,F2,IVECII(I),JLABI,JPARI, &
               -ENG, -AC, -OSCC, -SA, DSQRT(-SA)*DBLE(ISIGNA) !, ASFA
            WRITE (24, 301) -AB, -OSCB, -SB, DSQRT(-SB)*DBLE(ISIGNB) !, ASFB
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f1'
                  F2 = 'f2'
               END IF
               WRITE (24,300) F1,IVECII(I),JLABI,JPARI,F2,IVECFF(J), &
                     JLABJ,JPARJ,ENG,AC*IATJPOFF(J)/IATJPOII(I),     &
                     OSCC,SA,DSQRT(SA)*DBLE(ISIGNA)
               WRITE (24,301) AB*IATJPOFF(J)/IATJPOII(I),OSCB,SB,DSQRT(SB)*DBLE(ISIGNB)
            END IF
            LINES = LINES + 3
         ENDIF
!
      ELSE
!
!   Magnetic multipoles
!
!   In atomic units
!
         IF (ASFA < 0) THEN
            ISIGNA = -1
         ELSE
            ISIGNA = 1
         ENDIF
         AMS = ASFA**2
         AM = AMS*FAAU
         BM = AMS*FBAU
         OSCM = AMS*FOSC
         SA = OSCM*GGFACTOR*CVAC*CVAC*4
!GG-end
!
!   Convert to SI units if option 5 not set
!
         IF (.NOT.LTC(7)) THEN
            AM = AM*FASI
            BM = AM*FBSI
         ENDIF
!
!   Accumulate total of A coefficients
!
         TOTC(J) = TOTC(J) + AM
         TOTB(J) = TOTB(J) + AM
!
!   Print information if AM is greater than CUTOFF
!
         IF (ABS(AM) >= CUTOFF) THEN
            IF (-1.d+9 < ENG .and. ENG .LT. -1.0d-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f2'
                  F2 = 'f1'
               END IF
               WRITE (24,302) F1,IVECFF(J),JLABJ,JPARJ,F2,IVECII(I),&
                     JLABI,JPARI,-ENG,-AM,-OSCM,-SA,DSQRT(-SA)*DBLE(ISIGNA)
!     :               ,ASFA (M-value)
!            WRITE (24,*) 'Relative sign',ISIGNA
            ELSE IF ( 1.d+9 .GT. ENG .AND. ENG .GT. 1.D-2) THEN
               IF (LSAME) THEN
                  F1 = 'f '
                  F2 = 'f '
               ELSE
                  F1 = 'f1'
                  F2 = 'f2'
               END IF
               WRITE (24,302) F1,IVECII(I),JLABI,JPARI,F2,IVECFF(J),&
                     JLABJ,JPARJ,ENG,AM*IATJPOFF(J)/IATJPOII(I),    &
                     OSCM,SA,DSQRT(SA)*DBLE(ISIGNA)
!     :              ,ASFA (M-value)
            END IF
            LINES = LINES+2
         ENDIF
!
      ENDIF
!
      RETURN
!
!cjb  format for highly charged ions         F11.2 -> F13.2
  300 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F13.2,' C',1P,  &
           4D13.5)
  301 FORMAT(42X,' B',1P,4D13.5)
!cjb  format for highly charged ions         F11.2 -> F13.2
  302 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,0P,F13.2,' M',1P,  &
           4D13.5)
      RETURN
!
      END SUBROUTINE PRINTA
