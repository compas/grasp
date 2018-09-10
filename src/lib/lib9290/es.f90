!***********************************************************************
!                                                                      *
      SUBROUTINE ES(F, S2F, S3F) 
!                                                                      *
!   Evaluate the sum of the series                                     *
!                                                                      *
!                       infinity      n              k                 *
!              S  (F) =   Sum     (-1)  exp (n*F) / n                  *
!               k        n = 0                                         *
!                                                                      *
!   for k = 2, 3 to machine precision. This is a utility subroutine,   *
!   called by SUBROUTINEs NUCPOT and NCHARG.                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 28 Sep 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:25   2/14/04  
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
      REAL(DOUBLE), INTENT(IN) :: F 
      REAL(DOUBLE), INTENT(OUT) :: S2F 
      REAL(DOUBLE), INTENT(OUT) :: S3F 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: N 
      REAL(DOUBLE) :: FASE, EN, OBN, ENF, TERM2, TERM3, S2LAST 
!-----------------------------------------------
!
      N = 0 
      S2F = 0.0D00 
      S3F = 0.0D00 
      FASE = 1.0D00 
      N = N + 1 
      EN = DBLE(N) 
      OBN = 1.0D00/EN 
      FASE = -FASE 
      ENF = EXP(EN*F) 
      TERM2 = FASE*ENF*OBN*OBN 
      TERM3 = TERM2*OBN 
      S2LAST = S2F 
      S2F = S2F + TERM2 
      S3F = S3F + TERM3 
      DO WHILE(ABS(S2F) /= ABS(S2LAST)) 
         N = N + 1 
         EN = DBLE(N) 
         OBN = 1.0D00/EN 
         FASE = -FASE 
         ENF = EXP(EN*F) 
         TERM2 = FASE*ENF*OBN*OBN 
         TERM3 = TERM2*OBN 
         S2LAST = S2F 
         S2F = S2F + TERM2 
         S3F = S3F + TERM3 
      END DO 
!
      RETURN  
      END SUBROUTINE ES 
