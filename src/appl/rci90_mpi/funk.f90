!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION FUNK (X, N)
!                                                                      *
!   This  function  evaluates the KN(X) functions using the analytic   *
!   functions defined in tables 1 and 3 of Fullerton and Rinker.       *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      REAL(DOUBLE), INTENT(IN) :: X
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: NP
      INTEGER :: K, NN, I
      REAL(DOUBLE), DIMENSION(10,4) :: P
      REAL(DOUBLE), DIMENSION(2,4) :: B
      REAL(DOUBLE), DIMENSION(3,4) :: C
      REAL(DOUBLE), DIMENSION(5,4) :: D, E
      REAL(DOUBLE) :: XN, SUM, X2, BSUM, CSUM, DSUM, ESUM, XM
!-----------------------------------------------
!
!
      DATA P/ 8.8357293375D-01, -2.8259817381D-01, -5.8904879578D-01, &
         1.2500133434D-01, -3.2729913852D-02, 8.2888574511D-03, &
         -1.0327765800D-05, 6.3643668900D-05, 0.0D00, 0.0D00, -7.1740181754D-01&
         , 1.1780972274D00, -3.7499963087D-01, 1.3089675530D-01, &
         -3.8258286439D-02, -2.4297287300D-05, -3.5920148670D-04, &
         -1.7170090700D-05, 0.0D00, 0.0D00, 9.9999999987D-01, 1.9770200000D-08&
         , -7.5000050190D-01, 7.8540306316D-01, -3.4988601655D-01, &
         6.4596333000D-05, -9.8189080747D-03, 8.6513145800D-05, &
         -2.3969236620D-04, 0.0D00, 6.0000000002D00, -6.4305200000D-08, &
         2.1049413000D-06, -2.6711271500D-05, -1.3705236152D-01, &
         -6.3476104090D-04, -7.8739801501D-02, -1.9641740173D-03, &
         -3.4752369349D-03,  - 7.3145316220D-04/
!
      DATA B/ -3.19999594323D+02, 2.53900995981D00, -6.40514843293D+01, &
         7.11722714285D-01, 5.19010136460D+03, 8.28495496200D+01, &
         3.18150793824D+02, 4.33898867347D+01/
!
      DATA C/ -3.19999594333D+02, 2.53901020662D+00, 0.0D00, 6.40514843287D+01&
         , -7.11722686403D-01, 8.042207748D-04, 2.76805406060D+04, &
         -3.27039477790D+02, 0.0D00, 8.48402116837D+02, -2.56939867765D+01, &
         3.20844906346D-01/
!
      DATA D/ 5.018065179D+00, 7.151891262D+01, 2.116209929D+02, &
         3.140327478D+01, -1.0D00, 2.172386409D+02, 1.643364528D+03, &
         2.122244512D+03, -4.512004044D+01, 1.0D00, 8.540770444D00, &
         6.076242766D+01, 9.714630584D+01, 3.154973593D+01, 1.0D00, &
         5.9243015865D-01, 2.0596312871D00, 3.7785190424D00, 3.5614853214D00, &
         1.0D00/
!
      DATA E/ 2.669207401D+00, 5.172549669D+01, 2.969809720D+02, &
         5.364324164D+02, 1.535335924D+02, 1.155589983D+02, 1.292191441D+03, &
         3.831198012D+03, 2.904410075D+03, 0.0D00, 4.543392478D00, &
         3.514920169D+01, 6.019668656D+01, 8.468839579D00, 0.0D00, &
         3.1511867816D-01, 3.473245222D-01, 3.8791936870D-02, -1.3059741497D-03&
         , 0.0D00/
!
      DATA NP/ 8, 8, 9, 10/
!
      IF (X /= 0.0D00) THEN
         IF (N<0 .OR. N==2 .OR. N==4 .OR. N>5) GO TO 98
         IF (N - 3 <= 0) THEN
            IF (N - 3 /= 0) THEN
               K = N + 1
               XN = 1.0D00
               GO TO 4
            ENDIF
            K = N
            XN = 1.0D00/X**2
            GO TO 4
         ENDIF
         K = N - 1
         XN = 1.0D00/X**4
    4    CONTINUE
         IF (X <= 1.0D00) THEN
!
!   Calculate function for X < = 1
!
            NN = NP(K)
            SUM = 0.0D00
            DO I = 1, NN
               SUM = SUM + P(I,K)*XN
               XN = XN*X
            END DO
            X2 = X*X
            BSUM = B(1,K) + X2*(B(2,K)+X2)
            CSUM = C(1,K) + X2*(C(2,K)+X2*C(3,K))
            GO TO (6,8,7,8) K
    6       CONTINUE
            BSUM = BSUM*X
            GO TO 8
    7       CONTINUE
            BSUM = BSUM*X2
    8       CONTINUE
            SUM = SUM + BSUM*DLOG(X)/CSUM
            FUNK = SUM
            RETURN
         ENDIF
!
!   Calculate function for X > 1
!
         XN = 1.0D00
         DSUM = 0.0D00
         ESUM = 0.0D00
         DO I = 1, 5
            DSUM = DSUM + D(I,K)*XN
            ESUM = ESUM + E(I,K)*XN
            XN = XN/X
         END DO
         XM = -X
         SUM = DSUM*EXP(XM)/(ESUM*DSQRT(X**3))
         FUNK = SUM
         RETURN
      ENDIF
      IF (N /= 0) GO TO 99
      FUNK = P(1,1)
      RETURN
!
!   Error section
!
   98 CONTINUE
      WRITE (*, 302)
      STOP
   99 CONTINUE
      WRITE (*, 301)
      STOP
!
  301 FORMAT(/,' Attempt to calculate FUNK (0,N) for N > 0')
  302 FORMAT(/,' Attempt to calculate FUNK (X,N) for N other than',&
         ' 0, 1, 3 and 5.')
      RETURN
!
      END FUNCTION FUNK
