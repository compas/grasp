!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION FUNL (X, K)
!                                                                      *
!   This  function  evaluates the LK(X) functions using the analytic   *
!   functions defined  in table 5  and equations  (20) and  (21)  of   *
!   Fullerton and Rinker.                                              *
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
      INTEGER, INTENT(IN) :: K
      REAL(DOUBLE), INTENT(IN) :: X
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K1, I
      REAL(DOUBLE), DIMENSION(6,2) :: F
      REAL(DOUBLE), DIMENSION(3,2) :: G
      REAL(DOUBLE), DIMENSION(2,2) :: H
      REAL(DOUBLE) :: A, B, SUM, XN, X2, SUMG, SUMH, XM
!-----------------------------------------------
!
!
      DATA F/ 2.008188D00, -2.397605D00, 1.046471D00, -3.670660D-01, &
         6.374000D-02, -3.705800D-02, 1.646407D00, -2.092942D00, 9.623100D-01, &
         -2.549600D-01, 1.644040D-01, 0.0D00/
!
      DATA G/ 7.51198D-01, 1.38889D-01, 2.0886D-02, 1.37691D-01, -4.16667D-01, &
          - 9.7486D-02/
!
      DATA H/ -4.44444D-01, -3.472D-03, 4.44444D-01, 1.7361D-02/
!
      DATA A, B/ 2.2D00,  - 1.72D00/
!
      IF (K>=0 .AND. K<=1) THEN
         IF (X <= 2.0D00) THEN
            IF (X == 0.0D00) GO TO 6
!
!   Use rational approximation for X < 2
!
            K1 = K + 1
            SUM = 0.0D00
            XN = 1.0D00
            DO I = 1, 6
               SUM = SUM + XN*F(I,K1)
               XN = XN*X
            END DO
            X2 = X*X
            SUMG = G(1,K1) + X2*(G(2,K1)+X2*G(3,K1))
            SUMH = H(1,K1) + X2*X2*H(2,K1)
            XN = DLOG(X)
            SUMG = XN*(SUMG + XN*SUMH)
            IF (K /= 0) THEN
               SUM = SUM + SUMG
               GO TO 7
            ENDIF
            SUM = SUM + X*SUMG
            GO TO 7
         ENDIF
!
         SUM = A + B/X
         IF (K /= 0) SUM = SUM + (SUM + B/X)/X
         SUM = SUM/X
         XM = -X
         SUM = SUM*DEXP(XM)
         GO TO 7
    6    CONTINUE
         IF (K == 1) GO TO 98
         SUM = F(1,1)
    7    CONTINUE
         FUNL = SUM
         RETURN
!
!   Error section
!
   98    CONTINUE
         WRITE (*, 302)
         STOP
      ENDIF
      WRITE (*, 301)
      STOP
!
  301 FORMAT(/,'FUNL: K must be either 0 or 1')
  302 FORMAT(/,'FUNL: Attempt to calculate function for'/,&
         ' zero argument and K value of 1')
      RETURN
!
      END FUNCTION FUNL
