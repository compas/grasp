*
*     ------------------------------------------------------------------
*       B D A T A RYSYS
*     ------------------------------------------------------------------
*
      SUBROUTINE BDATA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C  from library libang.a
C      CALL INITT
C  from library libsqlsf1.a
      call factrl(32)
      CALL GLCONS
      CALL RIBLS
      CALL TERMLS
C  from library libsqlsf2.a
C	CALL TRMF
C	CALL TRMF15
      RETURN
      END
*
*     -------------------------------------------------------------
*      B L O C K   D A T A   T E R M L S
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
CGGV      BLOCK DATA TERMLS
      SUBROUTINE TERMLS
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MT/MT(40)
      COMMON/SKMT2/MTF(8),MT3(90)
      DATA MT/00100,10000,00300,20102,00104,30000,10202,10004,00500,
     *00100,20302,20102,40104,20104,00304,00104,20306,20106,00106,
     *20108,00308,00108,20110,00112,50000,10000,30202,10202,10404,
     *10204,30004,10004,30206,10206,10006,10208,30008,10008,10210,
     *10012/
      DATA MTF/60106,70000,50202,50004,50206,50008,50210,50012/
      DATA MT3/080108,
     *090000,070202,070004,070206,070008,070210,070012,070214,070016,
     *100110,
     *110000,090202,090004,090206,090008,090210,090012,090214,090016,
     *090218,090020,
     *120112,
     *130000,110202,110004,110206,110008,110210,110012,110214,110016,
     *110218,110020,110222,110024,
     *140114,
     *150000,130202,130004,130206,130008,130210,130012,130214,130016,
     *130218,130020,130222,130024,130226,130028,
     *160116,
     *170000,150202,150004,150206,150008,150210,150012,150214,150016,
     *150218,150020,150222,150024,150226,150028,150230,150032,
     *180118,
     *190000,170202,170004,170206,170008,170210,170012,170214,170016,
     *170218,170020,170222,170024,170226,170028,170230,170032,170234,
     *170036/
      END
*
*     -------------------------------------------------------------
*      S I X J
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | I/2  J/2  K/2 |                                            *
*     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville             October 1996   *
*
*
      SUBROUTINE SIXJ(I,J,K,L,M,N,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL SAVE
      DIMENSION RACA(0:4,0:4,0:4,0:4,0:4,0:4)
      DATA RACA/15625*1.D-20/
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DATA UNDEF/1.D-20/
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(I,J,K,L,M,N).EQ.0) RETURN
      ENDIF
      SAVE = .FALSE.
      IF (MAX0(I,J,K,L,M,N).LE.4) THEN
        SI = RACA(I,J,K,L,M,N)
        IF (SI .EQ. UNDEF) THEN
          SAVE = .TRUE.
        ELSE
          RETURN
        END IF
      END IF
      CALL GRACAH1(I,J,M,L,K,N,SI)
      IF(MOD(I+J+M+L,4).NE.0)SI=-SI
      IF (SAVE) RACA(I,J,K,L,M,N) =SI
      RETURN
      END
*
*     -------------------------------------------------------------
*      B L O C K   D A T A   R I B L S
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville             October 1996   *
*
CGGV      BLOCK DATA RIBLS
      SUBROUTINE RIBLS
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      COMMON/RIBOLSF/IMPTLSF(8),IMGTLSF(8),IMPNLSF(8),IMGNLSF(8)
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
      COMMON/RIBOLS3/IMPTLS3(90),IMGTLS3(90),IMPNLS3(90),IMGNLS3(90)
      DATA IMPTLS/1,2,3*3,3*6,16*9,16*25/
      DATA IMGTLS/1,2,3*5,3*8,16*24,16*40/
      DATA IMPNLS/2,1,3*6,3*3,16*25,16*9/
      DATA IMGNLS/2,1,3*8,3*5,16*40,16*24/
      DATA IMPTLSF/301,7*302/
      DATA IMGTLSF/301,7*308/
      DATA IMPNLSF/302,7*301/
      DATA IMGNLSF/308,7*301/
      DATA IMPTF/119*1,119*120/
      DATA IMGTF/119*119,119*238/
      DATA IMPNF/119*120,119*1/
      DATA IMGNF/119*238,119*119/
      DATA IMPTLS3/1,9*2,11,11*12,23,13*24,37,15*38,53,17*54,71,
     *19*72/
      DATA IMGTLS3/1,9*10,11,11*22,23,13*36,37,15*52,53,17*70,71,
     *19*90/
      DATA IMPNLS3/2,9*1,12,11*11,24,13*23,38,15*37,54,17*53,72,
     *19*71/
      DATA IMGNLS3/10,9*1,22,11*11,36,13*23,52,15*37,70,17*53,90,
     *19*71/
      END
*
*     -------------------------------------------------------------
*      N I N E L S
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2  J2/2  J3/2 |                                        *
*     |  L1/2  L2/2  L3/2 |                                        *
*     |  K1/2  K2/2  K3/2 |                                        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius, LITHUANIA                              January 1997 *
*
      SUBROUTINE NINELS(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      A=ZERO
      IF(I.EQ.1) THEN
        IN=0
        IF(ITTK(J1,J2,J3).EQ.0)RETURN
        IF(ITTK(L1,L2,L3).EQ.0)RETURN
        IF(ITTK(K1,K2,K3).EQ.0)RETURN
        IF(ITTK(J1,L1,K1).EQ.0)RETURN
        IF(ITTK(J2,L2,K2).EQ.0)RETURN
        IF(ITTK(J3,L3,K3).EQ.0)RETURN
        IN=1
        RETURN
      ENDIF
      IN=1
      CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
      RETURN
      END
*
*     -------------------------------------------------------------
*      N I N E
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2  J2/2  J3/2 |                                        *
*     |  L1/2  L2/2  L3/2 |                                        *
*     |  K1/2  K2/2  K3/2 |                                        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IF(I.EQ.1) THEN
        IN=0
        IF(ITTK(J1,J2,J3).EQ.0)RETURN
        IF(ITTK(L1,L2,L3).EQ.0)RETURN
        IF(ITTK(K1,K2,K3).EQ.0)RETURN
        IF(ITTK(J1,L1,K1).EQ.0)RETURN
        IF(ITTK(J2,L2,K2).EQ.0)RETURN
        IF(ITTK(J3,L3,K3).EQ.0)RETURN
        IN=1
        RETURN
      ENDIF
C
        N1=IABS(J1-K3)
        N2=IABS(L3-J2)
        N3=IABS(L1-K2)
        N4=IABS(J2-L3)
        N5=IABS(K2-L1)
        N6=IABS(J1-K3)
        MAX=MAX0(N1,N2,N3,N4,N5,N6)
        N1=J1+K3
        N2=L3+J2
        N3=J2+L3
        N4=K2+L1
        N5=J1+K3
        N6=L1+K2
        MIN=MIN0(N1,N2,N3,N4,N5,N6)
        IN=1
        AA=ZERO
        DO 1 IX=MAX,MIN,2
          CALL SIXJ(J1,J2,J3,L3,K3,IX,0,S1)
          CALL SIXJ(L1,L2,L3,J2,IX,K2,0,S2)
          CALL SIXJ(K1,K2,K3,IX,J1,L1,0,S3)
          X=S1*S2*S3*DBLE(IX+1)
          IF(MOD(IX,2).NE.0) X=-X
          AA=X+AA
    1   CONTINUE
      RETURN
      END
*
*     -----------------------------------------------------------------
*      G R A C A H 1
*     -----------------------------------------------------------------
*
      SUBROUTINE GRACAH1(I,J,K,L,M,N,RAC)
      IMPLICIT REAL*8(A-H,O-Z)
*
*      SUBROUTINE TO CALCULATE RACAH COEFFICIENTS.
*      THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE.
*      WRITTEN BY N. S. SCOTT
*      Modified by C. Froese Fischer, March 11, 1988 to use
*         table look-up
*      and
*                  G. Gaigalas, September 14, 1997
*
      COMMON/FACT/GAM(100)
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
*
      J1 = I+J+M
      J2 = K+L+M
      J3 = I+K+N
      J4 = J+L+N
      IF (MOD(J1,2) .EQ. 0  .AND.  MOD(J2,2) .EQ. 0   .AND.
     :    MOD(J3,2) .EQ. 0  .AND.  MOD(J4,2) .EQ. 0 )  THEN
          J1 = J1/2
          J2 = J2/2
          J3 = J3/2
          J4 = J4/2
          IF (MAX(I,J,M) .LE. J1 .AND.  MAX(K,L,M) .LE. J2  .AND.
     :        MAX(I,K,N) .LE. J3 .AND.  MAX(J,L,N) .LE. J4  )  THEN
              J5 = (I+J+K+L)/2
              J6 = (I+L+M+N)/2
              J7 = (J+K+M+N)/2
              NUMIN = MAX(J1, J2, J3, J4) + 1
              NUMAX = MIN(J5, J6, J7)     + 1
              RAC = ONE
              ICOUNT = 0
              DO 10 KK = NUMIN+1,NUMAX
                KI = NUMAX - ICOUNT
                RAC = ONE - (RAC*(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2)))/
     :                   ((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4))
                ICOUNT = ICOUNT+1
  10          CONTINUE
              RAC = RAC*EXP(
     :              (GAM(NUMIN+1) - GAM(NUMIN-J1) - GAM(NUMIN-J2) -
     :               GAM(NUMIN-J3) - GAM(NUMIN-J4) - GAM(J5+2-NUMIN)-
     :               GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN)) +
     :              (GAM(J1+1-I)+GAM(J1+1-J)+GAM(J1+1-M)-GAM(J1+2) +
     :               GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)-GAM(J2+2) +
     :               GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2) +
     :               GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))/TWO)
              IF (MOD(J5+NUMIN,2) .EQ. 0) RAC = -RAC
          ELSE
              RAC = ZERO
          END IF
      ELSE
         RAC = ZERO
      END IF
      RETURN
      END
*
*     -------------------------------------------------------------
*      I X J T I K
*     -------------------------------------------------------------
*
*                                                                  *
*     CHESKED TRIANGULAR CONDITIONS FOR 6j COEFFICIENT             *
*                                                                  *
*     | I/2  J/2  K/2 |            IXJTIK=1 - IF NOT SATISFY       *
*     | L/2  M/2  N/2 |            IXJTIK=0 - IN OVER CASES        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      FUNCTION IXJTIK(I,J,K,L,M,N)
      IXJTIK=0
      IF(ITTK(I,J,K).EQ.0)RETURN
      IF(ITTK(I,M,N).EQ.0)RETURN
      IF(ITTK(L,J,N).EQ.0)RETURN
      IF(ITTK(L,M,K).EQ.0)RETURN
      IXJTIK=1
      RETURN
      END
*
*     -------------------------------------------------------------
*      I T T K
*     -------------------------------------------------------------
*
*                                                                  *
*     CHESKED TRIANGULAR CONDITIONS FOR   I/2, J/2, K/2.           *
*     I+J>=K, I+K>=J, J+K>=I,                                      *
*     I/2+J/2+K/2 - WHOLE NUMBER                                   *
*     ITTK=1 -   IF NOT SATISFY                                    *
*     ITTK=0 -   IN OVER CASES                                     *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      FUNCTION ITTK(I,J,K)
      ITTK=0
      IF(IABS(I-J).GT.K)RETURN
      IF(I+J.LT.K)RETURN
      IF(MOD(I+J+K,2).NE.0)RETURN
      ITTK=1
      RETURN
      END
*
*     ---------------------------------------------------------------
*      I T R E X G
*     ---------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville             October 1996   *
*
      FUNCTION ITREXG(I1,I2,I3,I4,K)
      J=MAX0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG=J
      RETURN
      END
*
*     ------------------------------------------------------------------
*       B L O C K   D A T A    G L C O N S
*     ------------------------------------------------------------------
*
CGGV      BLOCK DATA GLCONS
      SUBROUTINE GLCONS
C
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
C     SET GLOBAL REAL CONSTANTS
C
      DATA ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS/
     :     0.0D 00,0.1D 00,0.5D 00,
     :     1.0D 00,2.0D 00,3.0D 00,
     :     4.0D 00,
     :     7.0D 00,1.1D 01,1.0D-08/
C
      END
*
*     -----------------------------------------------------------------
*           F A C T R L
*     -----------------------------------------------------------------
*
*
      SUBROUTINE FACTRL(NFACT)
*
*      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON/FACT/GAM(100)
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
*
      GAMMA=ONE
      GAM(1) = ZERO
      DO 1 I=1,NFACT-1
         GAMMA=I*GAMMA
         GAM(I+1) = DLOG(GAMMA)
    1 CONTINUE
      DO 20 I = NFACT+1,(100)
         X = I-1
         GAM(I) = GAM(I-1) + DLOG(X)
   20 CONTINUE
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE JJPER(IL1,IL2,IS1,IS2,IL,IS,IJ,IJ1,IJ2,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      CALL NINELS(IL1,IS1,IJ1,IL2,IS2,IJ2,IL,IS,IJ,1,INN,SN)
      IF(INN.NE.0) THEN
        CALL NINELS(IL1,IS1,IJ1,IL2,IS2,IJ2,IL,IS,IJ,0,INN,SN)
        REZ=DSQRT(DBLE((IL+1)*(IS+1)*(IJ1+1)*(IJ2+1)))*SN
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE JKPER(IL1,IL2,IS1,IS2,IL,IS,IJ,IJ1,IK,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ILES=IXJTIK(IL1,IL2,IL,IK,IS1,IJ1)
      IF(ILES.EQ.1) THEN
        ISES=IXJTIK(IS1,IS2,IS,IJ,IL,IK)
        IF(ISES.EQ.1) THEN
          CALL SIXJ(IL1,IL2,IL,IK,IS1,IJ1,0,SNL)
          CALL SIXJ(IS1,IS2,IS,IJ,IL,IK,0,SNS)
          REZ=DSQRT(DBLE((IJ1+1)*(IK+1)*(IL+1)*(IS+1)))*SNL*SNS
          IF(MOD(IL2-IS2+IJ1-IJ,4).NE.0) REZ=-REZ
        ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE LKPER(IL1,IL2,IS1,IS2,IL,IS,IJ,ILS,IK,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      IF(IL.NE.ILS) RETURN
      ISES=IXJTIK(IS1,IS2,IS,IJ,IL,IK)
      IF(ISES.NE.0) THEN
         CALL SIXJ(IS1,IS2,IS,IJ,IL,IK,0,SN)
         REZ=DSQRT(DBLE((IK+1)*(IS+1)))*SN
         IF(MOD(IL+IS1+IS2+IJ,4).NE.0) REZ=-REZ
      ENDIF
      RETURN
      END

*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
c	SUBROUTINE JJJKPER(I1,I2,IK,ISI,IJ,IJI,IJIM1,ILI,REZ)
      SUBROUTINE JJJKPER(IK,ISI,IJ,IJI,IJIM1,ILI,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
c      IF(I1.NE.I2) RETURN
      ISES=IXJTIK(IK,ISI,IJ,IJI,IJIM1,ILI)
      IF(ISES.NE.0) THEN
         CALL SIXJ(IK,ISI,IJ,IJI,IJIM1,ILI,0,SN)
         REZ=DSQRT(DBLE((IK+1)*(IJI+1)))*SN
         IF(MOD(IJIM1+ILI+ISI+IJ,4).NE.0) REZ=-REZ
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE LS3PER(L_12, L_3, L_123, L_4, L_34, L,
     :                             IS_12,IS_3,IS_123,IS_4,IS_34,IS,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(L_12,L_3,L_123,L_4,L,L_34)
      IF(ISES.NE.0) THEN
         ISES=IXJTIK(IS_12,IS_3,IS_123,IS_4,IS,IS_34)
         IF(ISES.NE.0) THEN
            CALL SIXJ(L_12,L_3,L_123,L_4,L,L_34,0,SNL)
            CALL SIXJ(IS_12,IS_3,IS_123,IS_4,IS,IS_34,0,SNS)
            REZ=DSQRT(DBLE((L_123+1)*(L_34+1)*(IS_123+1)*(IS_34+1)))
            REZ = REZ*SNL*SNS
            IF(MOD(L_3+L_4+L_12+L+IS_3+IS_4+IS_12+IS,4).NE.0) REZ=-REZ
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE LSJ3PER(L_12, L_3, L_123, L_4, L_34, L,
     :                  IS_12,IS_3,IS_123,IS_4,IS_34,IS,J_12,J_34,J,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(L_12,L_3,L_123,L_4,L,L_34)
      IF(ISES.NE.0) THEN
         ISES=IXJTIK(IS_12,IS_3,IS_123,IS_4,IS,IS_34)
         IF(ISES.NE.0) THEN
            CALL NINELS(L_12,L_34,L,IS_12,IS_34,IS,J_12,J_34,J,1,INN,SN)
            IF(INN.NE.0) THEN
               CALL NINELS
     :                 (L_12,L_34,L,IS_12,IS_34,IS,J_12,J_34,J,0,INN,SN)
               CALL SIXJ(L_12,L_3,L_123,L_4,L,L_34,0,SNL)
               CALL SIXJ(IS_12,IS_3,IS_123,IS_4,IS,IS_34,0,SNS)
               REZ=DSQRT(DBLE((L_123+1)*(L_34+1)*(IS_123+1)*(IS_34+1)*
     :                        (L+1)*(IS+1)*(J_12+1)*(J_34+1)))
               REZ = REZ*SNL*SNS*SN
               IF(MOD(L_3+L_4+L_12+L+IS_3+IS_4+IS_12+IS,4).NE.0)REZ=-REZ
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE LK3PER(L_12, L_3, L_123, L_4, L_34, L,
     :                  IS_12,IS_3,IS_123,IS_4,IS_34,IS,K,J,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(L_12,L_3,L_123,L_4,L,L_34)
      IF(ISES.NE.0) THEN
         ISES=IXJTIK(IS_12,IS_3,IS_123,IS_4,IS,IS_34)
         IF(ISES.NE.0) THEN
            ISES=IXJTIK(IS_34,IS_12,IS,L,J,K)
            IF(ISES.NE.0) THEN
               CALL SIXJ(IS_34,IS_12,IS,L,J,K,0,SN)
               CALL SIXJ(L_12,L_3,L_123,L_4,L,L_34,0,SNL)
               CALL SIXJ(IS_12,IS_3,IS_123,IS_4,IS,IS_34,0,SNS)
               REZ=DSQRT(DBLE((L_123+1)*(L_34+1)*(IS_123+1)*(IS_34+1)*
     :                       (K+1)*(IS+1)))
               REZ = REZ*SNL*SNS*SN
               IF(MOD(
     :            L_3+L_4+L_12+L+IS_3+IS_4+IS_12+IS+IS_12+IS_34+L+J,4)
     :                                                  .NE.0) REZ=-REZ
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE JK3PER(L_12, L_3, L_123, L_4, L_34, L,
     :                  IS_12,IS_3,IS_123,IS_4,IS_34,IS,J_12,K,J,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(L_12,L_3,L_123,L_4,L,L_34)
      IF(ISES.NE.0) THEN
         ISES=IXJTIK(IS_12,IS_3,IS_123,IS_4,IS,IS_34)
         IF(ISES.NE.0) THEN
            ISES=IXJTIK(IS_12,L_12,J_12,L_34,K,L)
            IF(ISES.NE.0) THEN
               ISES=IXJTIK(IS_12,IS_34,IS,J,L,K)
               IF(ISES.NE.0) THEN
                 CALL SIXJ(IS_12,L_12,J_12,L_34,K,L,0,SN1)
                 CALL SIXJ(IS_12,IS_34,IS,J,L,K,0,SN2)
                 CALL SIXJ(L_12,L_3,L_123,L_4,L,L_34,0,SNL)
                 CALL SIXJ(IS_12,IS_3,IS_123,IS_4,IS,IS_34,0,SNS)
                 REZ=DSQRT(DBLE((L_123+1)*(L_34+1)*(IS_123+1)*(IS_34+1)*
     :                      (J_12+1)*(L+1))*(IS+1)*(K+1))
                 REZ = REZ*SNL*SNS*SN1*SN2
                 IF(MOD(
     :           L_3+L_4+L_12+L+IS_3+IS_4+IS_12+IS+
     :           J_12+L_34+K+IS_34+K+J,4)
     :                                                 .NE.0) REZ=-REZ
              ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE cLSJ3PER(L_1,L_12, L_2, L_3, L_23, L,
     :                  IS_1,IS_12,IS_2,IS_3,IS_23,IS,J_1,J_12,J,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(L_1,L_2,L_12,L_3,L,L_23)
      IF(ISES.NE.0) THEN
         ISES=IXJTIK(IS_1,IS_2,IS_12,IS_3,IS,IS_23)
         IF(ISES.NE.0) THEN
            CALL NINELS(L_1,IS_1,J_1,L_23,IS_23,J_12,L,IS,J,1,INN,SN)
            IF(INN.NE.0) THEN
               CALL NINELS
     :                 (L_1,IS_1,J_1,L_23,IS_23,J_12,L,IS,J,0,INN,SN)
               CALL SIXJ(L_1,L_2,L_12,L_3,L,L_23,0,SNL)
               CALL SIXJ(IS_1,IS_2,IS_12,IS_3,IS,IS_23,0,SNS)
               REZ=DSQRT(DBLE((L_12+1)*(L_23+1)*(IS_12+1)*(IS_23+1)*
     :                        (L+1)*(IS+1)*(J_1+1)*(J_12+1)))
               REZ = REZ*SNL*SNS*SN
               IF(MOD(L_2+L_3+L_1+L+IS_2+IS_3+IS_1+IS,4).NE.0)REZ=-REZ
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE LScjjPER(L12,IS12,L3,IS3,L,IS,J_1,J_2,J,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      CALL NINELS(IS3,L3,J_2,IS12,L12,J_1,IS,L,J,1,INN,SN)
      IF(INN.NE.0) THEN
         CALL NINELS(IS3,L3,J_2,IS12,L12,J_1,IS,L,J,0,INN,SN)
         REZ=DSQRT(DBLE((J_1+1)*(J_2+1)*(L+1)*(IS+1)))
         REZ = REZ*SN
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE jj2PER(L1,IS1,L2,IS2,L,IS,J_1,J_2,J,JM2,JP2,J_12,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(JM2,JP2,J_2,J,J_1,J_12)
      IF(ISES.NE.0) THEN
         CALL NINELS(IS2,L2,J_2,IS1,L1,J_1,IS,L,J,1,INN,SN)
         IF(INN.NE.0) THEN
            CALL NINELS(IS2,L2,J_2,IS1,L1,J_1,IS,L,J,0,INN,SN)
            CALL SIXJ(JM2,JP2,J_2,J,J_1,J_12,0,SNS)
            REZ=DSQRT(DBLE((J_1+1)*(J_2+1)*(L+1)*(IS+1)*
     :                             (J_2+1)*(J_12+1)))
            REZ = REZ*SN*SNS
            IF(MOD(JM2+JP2+J_1+J,4).NE.0)REZ=-REZ
         ENDIF
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*
*     -----------------------------------------------------------------
*
      SUBROUTINE jj3PER(L1,IS1,L2,IS2,L_12,IS_12,L3,IS3,L,IS,
     :           J_1,J_2,J_3,J,JM2,JP2,J_12,J_12S,JM3,JP3,J_123S,REZ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
C
      REZ=ZERO
      ISES=IXJTIK(JM2,JP2,J_2,J_12,J_1,J_12S)
      IF(ISES.EQ.0) RETURN
      ISES=IXJTIK(JM3,JP3,J_3,J,J_12,J_123S)
      IF(ISES.EQ.0) RETURN
      CALL NINELS(L1,IS1,J_1,L2,IS2,J_2,L_12,IS_12,J_12,1,INN,SN)
      IF(INN.EQ.0) RETURN
      CALL NINELS(L_12,L3,L,IS_12,IS3,IS,J_12,J_3,J,1,INN,SN)
      IF(INN.EQ.0) RETURN
      CALL NINELS(L1,IS1,J_1,L2,IS2,J_2,L_12,IS_12,J_12,0,INN,SN1)
      CALL NINELS(L_12,L3,L,IS_12,IS3,IS,J_12,J_3,J,0,INN,SN2)
      CALL SIXJ(JM2,JP2,J_2,J_12,J_1,  J_12S,0,SNS1)
      CALL SIXJ(JM3,JP3,J_3,J,   J_12,J_123S,0,SNS2)
      REZ = SN1*SN2*SNS1*SNS2
      REZ=REZ*DSQRT(DBLE((L_12+1)*(L+1)*(IS_12+1)*(IS+1)))
      REZ=REZ*DSQRT(DBLE((J_1+1)*(J_2+1)*(J_3+1)*(J_12+1)))
      REZ=REZ*DSQRT(DBLE((J_2+1)*(J_12S+1)*(J_3+1)*(J_123S+1)))

      IF(MOD(JM2+JP2+J_1+J_12+JM3+JP3+J_12+J,4).NE.0)REZ=-REZ
      RETURN
      END
*
*     -----------------------------------------------------------------
*       L V A L
*     -----------------------------------------------------------------
*
*     Modified by Gediminas Gaigalas,                September 1997
*
      INTEGER FUNCTION LVAL(SYMBOL)
      CHARACTER*1 SYMBOL
      CHARACTER*26 SET
      DATA SET/'spdfghiklmnoqSPDFGHIKLMNOQ'/
*
      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE .LE. 13) THEN
            LVAL = LOCATE - 1
         ELSE
            LVAL = LOCATE - 14
      ENDIF
      RETURN
      END
*
*     -----------------------------------------------------------------
*       C V A L
*     -----------------------------------------------------------------
*
*     Writen by Gediminas Gaigalas,                September 1999
*
*      SUBROUTINE CVAL(I,IVALUE,SYMBOL)
      CHARACTER*1 FUNCTION CVAL(I,IVALUE)
      IMPLICIT REAL*8(A-H,O-Z)
*      CHARACTER*1 SYMBOL
      CHARACTER*26 SET
      DATA SET/'spdfghiklmnoqSPDFGHIKLMNOQ'/
*
      IF (I .EQ. 1) THEN
         II = IVALUE+1
      ELSE IF (I .EQ. 2) THEN
         II = IVALUE/2
         II = II+14
      ELSE
        PRINT*, "error in CVAL"
        STOP
      ENDIF
*      SYMBOL = SET(II:II)
      CVAL = SET(II:II)
      RETURN
      END
*
*     -----------------------------------------------------------------
*       J V A L
*     -----------------------------------------------------------------
*
*     Writen by Gediminas Gaigalas,                September 1999
*
*      SUBROUTINE JVAL(IVALUE,SYMBOL)
      CHARACTER*4 FUNCTION JVAL(IVALUE)
      IMPLICIT REAL*8(A-H,O-Z)
*      CHARACTER*1 SYMBOL
      CHARACTER*64 SET1, SET2
      CHARACTER*8  SET3
      DATA SET1/
     :'1/2 1   3/2 2   5/2 3   7/2 4   9/2 5   11/26   13/27   15/28   '
     :/
      DATA SET2/
     :'17/29   19/210  21/211  23/212  25/213  27/214  29/215  31/216  '
     :/
      DATA SET3/'33/217  '/
*
      IF( IVALUE .EQ. 0) THEN
         JVAL = '0   '
         I = 1
      ELSE IF (IVALUE .LE. 16) THEN
         I = (IVALUE)*4
         JVAL = SET1(I-3:I)
      ELSE IF (IVALUE .LE. 32) THEN
         I = (IVALUE-16)*4
         JVAL = SET2(I-3:I)
      ELSE IF (IVALUE .LE. 34) THEN
         I = (IVALUE-32)*4
         JVAL = SET3(I-3:I)
      ELSE
         print*,"IVALUE =",IVALUE
         STOP ' ERROR in JVAL'
      END IF
      RETURN
      END
*
*     -----------------------------------------------------------------
*       L S J V A L
*     -----------------------------------------------------------------
*
*     Writen by Gediminas Gaigalas,                September 1999
*
*      SUBROUTINE JVAL(IVALUE,SYMBOL)
      CHARACTER*2 FUNCTION LSJVAL(IL,IJ)
      IMPLICIT REAL*8(A-H,O-Z)
*      CHARACTER*1 SYMBOL
      CHARACTER*52 SET
      DATA SET/
     :'s-s+p-p+d-d+f-f+g-g+h-h+i-i+k-k+l-l+m-m+n-n+o-o+q-q+'
     :/
*
      IF (IL .LE. 13) THEN
         I = IJ-2*IL
         IF(I .EQ. -1) THEN
           I = 4*IL + 1
         ELSE IF (I .EQ. 1) THEN
           I = 4*IL + 3
         ELSE IF (IL .EQ. 0 .AND. IJ .EQ. 0) THEN
           I = 1
         ELSE
          print*,"IL",IL,"  IJ =",IJ
          STOP ' ERROR in LSJVAL 1'
         END IF
         LSJVAL = SET(I:I+2)
      ELSE
         STOP ' ERROR in LSJVAL 2'
      END IF
      RETURN
      END
