!
!***********************************************************************
!                                                                      *
      SUBROUTINE PACKLS(M, EL, Q, COUPLE, STR) 
!                                                                      *
!   Subroutine written by Bin LIU                                      *
!                                                                      *
!   Rules for encoding                                                 *
!   1. All blanks deleted                                              *
!   2. If Qi=1, omit Qi                                                *
!   3. If Qi=1 or Qi>=4l+1, omit ALFAi                                 *
!   4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.          *
!                                                                      *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lval_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)            :: M 
!GG      CHARACTER(LEN=64), INTENT(OUT) :: STR 
      CHARACTER(LEN=164), INTENT(OUT) :: STR 
      INTEGER, INTENT(IN)            :: Q(*) 
      CHARACTER(LEN=3), INTENT(IN)   :: EL(*), COUPLE(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: FULL, CONST, J, I, K, N 
      CHARACTER(LEN=1) :: CH1 
      CHARACTER(LEN=3) :: CH3 
!-----------------------------------------------
!  FULL   :  4l+2
!  CONST  :   constant for converting lowercase to uppercase
!  CH*    :  temporary variables
!
      CONST = ICHAR('a') - ICHAR('A') 
      STR = ' ' 
      J = 1 
!
!   -----  begin to encode  -----
!
      DO I = 1, M 
         IF (Q(I) == 0) CYCLE  
         IF (EL(I)(1:1) == ' ') THEN 
            STR(J:J+1) = EL(I)(2:3) 
            K = 2 
         ELSE 
            STR(J:J+2) = EL(I) 
            K = 3 
            IF (EL(I)(3:3) == ' ') K = 2 
         ENDIF 
         CH1 = STR(J+1:J+1) 
         IF (CH1>='A' .AND. CH1<='Z') STR(J+1:J+1) = CHAR(ICHAR(CH1) + CONST) 
         FULL = 4*LVAL(CH1) + 2 
!
!  -----  convert Qi into character  -----
!
         J = J + K 
         N = Q(I) 
         IF (N > 14) STOP 'Too many electrons' 
!
!  -----  If Qi<>1, add Qi
!           If Qi<4l+1, add TERMi for the shell -----
!
         IF (N /= 1) THEN 
            STR(J:J) = '(' 
            J = J + 1 
            IF (N > 9) THEN 
               STR(J:J) = '1' 
               J = J + 1 
               N = N - 10 
            ENDIF 
            STR(J:J+1) = CHAR(ICHAR('0') + N)//')' 
            J = J + 2 
            IF (N<FULL - 1 .AND. M/=1) THEN 
               CH3 = COUPLE(I) 
               CH1 = CH3(2:2) 
               IF (CH1>='a' .AND. CH1<='z') CH3(2:2) = CHAR(ICHAR(CH1) - CONST) 
               STR(J:J+2) = CH3 
               J = J + 3 
            ENDIF 
         ENDIF 
!
!  -----  If i=1 or Qi=4l+2 and i<>m,
!           insert '.'; else _RESULTANTi.  -----
!
         IF (I/=1 .AND. N/=FULL .OR. M==I) THEN 
            CH3 = COUPLE(M+I-1) 
            CH1 = CH3(2:2) 
            IF (CH1>='a' .AND. CH1<='z') CH3(2:2) = CHAR(ICHAR(CH1) - CONST) 
            K = 2 
            IF (M == 1) K = 3 
            STR(J:J+K) = '_'//CH3(1:K) 
            J = J + K + 1 
         ENDIF 
         IF (I == M) CYCLE  
         STR(J:J) = '.' 
         J = J + 1 
      END DO 
!
!>>>>>    Because of a compiler error on the SUN, the following is
!         needed to have the string printed correctly.
!         STR = STR(1:J)
      RETURN  
      END SUBROUTINE PACKLS 
