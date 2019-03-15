!***********************************************************************
!                                                                      *
      SUBROUTINE NUCPOT
!                                                                      *
!   Evaluate the nuclear potential for point and Fermi models.         *
!                                                                      *
!   Call(s) to: [LIB92] ES.                                            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:02   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE DEBUG_C
      USE DEF_C,           ONLY: Z, PI
      USE GRID_C
      USE NPAR_C
      USE NPOT_C,          ONLY: ZZ, NNUC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE es_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NB3, NROWS, II, II1, II2, II3
      REAL(DOUBLE) :: C, A, ABC, TABC, ABC2, THABC2, ABC3, CBA, PI2, HPIAC2, &
         H3PHP, S2MCBA, S3MCBA, SABC3, DMSAS, EN, ZBN, RI, RMC, RMCBA, RBC, &
         S2RCBA, S3RCBA
      LOGICAL :: SET
!-----------------------------------------------
!
!
!   Point nucleus
!
      IF (NPARM == 0) THEN
!
         ZZ(:N) = Z
!
!   Fermi distribution
!
      ELSE IF (NPARM == 2) THEN
!
         C = PARM(1)
         A = PARM(2)
         ABC = A/C
         TABC = 2.0D00*ABC
         ABC2 = ABC*ABC
         THABC2 = 3.0D00*ABC2
         ABC3 = ABC2*ABC
         CBA = C/A
         PI2 = PI*PI
         HPIAC2 = 0.5D00*PI2*ABC2
         H3PHP = 1.5D00 + HPIAC2
         CALL ES ((-CBA), S2MCBA, S3MCBA)
         SABC3 = 6.0D00*ABC3
         DMSAS = -SABC3*S3MCBA
         EN = 1.0D00 + ABC2*PI2 + DMSAS
         ZBN = Z/EN
!
         SET = .FALSE.
         DO I = 1, N
            RI = R(I)
            RMC = RI - C
            RMCBA = RMC/A
            RBC = RI/C
            IF (RBC <= 1.0D00) THEN
               CALL ES (RMCBA, S2RCBA, S3RCBA)
               ZZ(I) = ZBN*(DMSAS + SABC3*S3RCBA + RBC*(H3PHP - THABC2*S2RCBA&
                   - 0.5D00*RBC*RBC))
            ELSE
               IF (.NOT.SET) THEN
                  NNUC = I
                  SET = .TRUE.
               ENDIF
               CALL ES ((-RMCBA), S2RCBA, S3RCBA)
               ZZ(I) = Z*(1.0D00 + THABC2*(RBC*S2RCBA + TABC*S3RCBA)/EN)
            ENDIF
         END DO
      ENDIF
!
      IF (LDBPR(2)) THEN
         WRITE (99, 300)
         NB3 = N/3
         IF (3*NB3 == N) THEN
            NROWS = NB3
         ELSE
            NROWS = NB3 + 1
         ENDIF
         DO II = 1, NROWS
            II1 = II
            II2 = II1 + NROWS
            II3 = II2 + NROWS
            IF (II3 <= N) THEN
               WRITE (99, 301) R(II1), ZZ(II1), R(II2), ZZ(II2), R(II3), ZZ(II3&
                  )
            ELSE IF (II2 <= N) THEN
               WRITE (99, 301) R(II1), ZZ(II1), R(II2), ZZ(II2)
            ELSE
               WRITE (99, 301) R(II1), ZZ(II1)
            ENDIF
         END DO
      ENDIF
!
      RETURN
!
  300 FORMAT(/,'From SUBROUTINE NUCPOT:'/,3(&
         ' -------- r -------- ----- -r*V(r) -----'))
  301 FORMAT(1P,6(1X,1D19.12))
      RETURN
!
      END SUBROUTINE NUCPOT
