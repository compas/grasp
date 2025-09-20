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
      INTEGER          :: FULL1, FULL, CONST, J, I, K, N, N1, NN, ICORE
      INTEGER          :: L
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
      ICORE = 0
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
         IF(I /= 1) FULL1 = FULL
         IF(I /= 1) N1 = N
         L = LVAL(CH1)
         FULL = 4*L + 2
!
!  -----  convert Qi into character  -----
!
         J = J + K
         N = Q(I)
         IF (N > 14) STOP 'Too many electrons'
         IF(I /= 1 .AND. FULL1 /= N1) ICORE = ICORE + FULL1 - N1
!
!  -----  If Qi<>1, add Qi
!           If Qi<4l+1, add TERMi for the shell -----
!
         IF (N /= 1) THEN
            STR(J:J) = '('
            J = J + 1
            NN = N
            IF (N > 9) THEN
               STR(J:J) = '1'
               J = J + 1
               NN = N - 10
            ENDIF
            STR(J:J+1) = CHAR(ICHAR('0') + NN)//')'
            J = J + 2
            IF(I/=1 .AND. ICORE<2 .AND. Q(I-1)<=1 .AND. N==2) THEN
               CONTINUE 
            ELSE IF(I/=1 .AND. ICORE==0 .AND. Q(I-1)==FULL1 .AND. N<3  &
                                                      .AND. M==I) THEN
               CONTINUE 
            ELSE IF(I/=1 .AND. ICORE==0 .AND. Q(I-1)==FULL1 .AND. L<2  &
                                                      .AND. M==I) THEN
               CONTINUE 
            ELSE IF(I/=1 .AND. ICORE==0 .AND. Q(I-1)==FULL1 .AND. L==2 &
                                          .AND. N > 7 .AND. M==I) THEN
               CONTINUE 
            ELSE IF(I/=1 .AND. ICORE==0 .AND. Q(I-1)==FULL1 .AND. L==3 &
                                         .AND. N > 11 .AND. M==I) THEN
               CONTINUE 
            ELSE IF (N<(FULL-1) .AND. M/=1) THEN
               CH3 = COUPLE(I)
               CH1 = CH3(2:2)
               IF (CH1>='a' .AND. CH1<='z') CH3(2:2) = CHAR(ICHAR(CH1) - CONST)
               IF(Len_Trim(CH3) == 2) THEN
                  STR(J:J+1) = TRIM(CH3)
                  J = J + 2
               ELSE IF(L == 2 .AND. N < 3) THEN
                  STR(J:J+1) = CH3(1:2)
                  J = J + 2
               ELSE IF(L == 2 .AND. N > 7) THEN
                  STR(J:J+1) = CH3(1:2)
                  J = J + 2
               ELSE IF(L == 2 .AND. N == 3) THEN
                  IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 2 .AND. N == 7) THEN
                  IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 2 .AND. N == 4) THEN
                  IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'S') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 2 .AND. N == 6) THEN
                  IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'S') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 2 .AND. N == 5) THEN
                  IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N < 3) THEN
                  STR(J:J+1) = CH3(1:2)
                  J = J + 2
               ELSE IF(L == 3 .AND. N > 11) THEN
                  STR(J:J+1) = CH3(1:2)
                  J = J + 2
               ELSE IF(L == 3 .AND. N == 3) THEN
                  IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N == 11) THEN
                  IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N == 4) THEN
                  IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'L') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'S') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N == 10) THEN
                  IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '3' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'L') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '1' .AND. CH3(2:2) == 'S') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N == 5) THEN
                  IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'M') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'L') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE IF(L == 3 .AND. N == 9) THEN
                  IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '4' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'M') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'L') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'K') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'I') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'H') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'G') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'F') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'D') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE IF(CH3(1:1) == '2' .AND. CH3(2:2) == 'P') THEN
                    STR(J:J+2) = CH3
                    J = J + 3
                  ELSE
                    STR(J:J+1) = CH3(1:2)
                    J = J + 2
                  END IF
               ELSE
                  STR(J:J+2) = CH3
                  J = J + 3
               END IF
            ENDIF
         ENDIF
!
!  -----  If i=1 or Qi=4l+2 and i<>m,
!           insert '.'; else _RESULTANTi.  -----
!
         IF(I/=1 .AND. Q(I-1)==FULL1 .AND. N==1 .AND. M/=I .AND.      &
                                                       ICORE<2) THEN
            CONTINUE 
         ELSE IF(I/=1 .AND. Q(I-1)==FULL1 .AND. N==FULL-1 .AND.       &
                                            M/=I .AND. ICORE<2) THEN
            CONTINUE 
         ELSE IF(I/=1 .AND. Q(I-1)==FULL1 .AND. ICORE==0 .AND.        &
                                                          M/=I) THEN
            CONTINUE 
         ELSE IF (I/=1 .AND. N/=FULL .OR. M==I) THEN
            CH3 = COUPLE(M+I-1)
            CH1 = CH3(2:2)
            IF (CH1>='a' .AND. CH1<='z') CH3(2:2)=CHAR(ICHAR(CH1)-CONST)
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
