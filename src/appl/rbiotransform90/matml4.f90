!
!     ------------------------------------------------------------------
!     M A T M L 4
!     ------------------------------------------------------------------
!
      SUBROUTINE MATML4(C, A, B, NCROW, NCCOL, NAROW, NACOL, NBROW, NBCOL, &
         ITRNSP) 
!
! MULTIPLY A AND B TO GIVE C
!
!     C = A * B             FOR ITRNSP = 0
!
!     C = A(TRANSPOSED) * B FOR ITRNSP = 1
!
!     C = A * B(TRANSPOSED) FOR ITRNSP = 2
!
!... JEPPE OLSEN, LAST REVISION JULY 24 1987
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE wrtmat_I 
      USE setvec_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCROW 
      INTEGER  :: NCCOL 
      INTEGER  :: NAROW 
      INTEGER  :: NACOL 
      INTEGER  :: NBROW 
      INTEGER  :: NBCOL 
      INTEGER , INTENT(IN) :: ITRNSP 
      REAL(DOUBLE)  :: C(NCROW,NCCOL) 
      REAL(DOUBLE)  :: A(NAROW,NACOL) 
      REAL(DOUBLE)  :: B(NBROW,NBCOL) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTEST, J, K, I 
      REAL(DOUBLE) :: BKJ, BJK 
!-----------------------------------------------
!
      NTEST = 0 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) 
         WRITE (6, *) ' A AND B MATRIX FROM MATML4 ' 
         WRITE (6, *) 
         CALL WRTMAT (A, NAROW, NACOL, NAROW, NACOL) 
         CALL WRTMAT (B, NBROW, NBCOL, NBROW, NBCOL) 
         WRITE (6, *) ' NCROW NCCOL NAROW NACOL NBROW NBCOL ' 
         WRITE (6, '(6I6)') NCROW, NCCOL, NAROW, NACOL, NBROW, NBCOL 
      ENDIF 
!
      CALL SETVEC (C, 0.0D0, NCROW*NCCOL) 
!
      IF (ITRNSP == 0) THEN 
         DO J = 1, NCCOL 
            DO K = 1, NBROW 
               BKJ = B(K,J) 
               C(:NCROW,J) = C(:NCROW,J) + A(:NCROW,K)*BKJ 
            END DO 
         END DO 
      ENDIF 
!
!
      IF (ITRNSP == 1) THEN 
!... C = A(T) * B
         DO J = 1, NCCOL 
            DO K = 1, NBROW 
               BKJ = B(K,J) 
               C(:NCROW,J) = C(:NCROW,J) + A(K,:NCROW)*BKJ 
            END DO 
         END DO 
      ENDIF 
!
      IF (ITRNSP == 2) THEN 
!... C = A*B(T)
         DO J = 1, NCCOL 
            DO K = 1, NBCOL 
               BJK = B(J,K) 
               C(:NCROW,J) = C(:NCROW,J) + A(:NCROW,K)*BJK 
            END DO 
         END DO 
      ENDIF 
!
!
      IF (NTEST /= 0) THEN 
         WRITE (6, *) 
         WRITE (6, *) ' C MATRIX FROM MATML4 ' 
         WRITE (6, *) 
         CALL WRTMAT (C, NCROW, NCCOL, NCROW, NCCOL) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE MATML4 
