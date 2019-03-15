!***********************************************************************
!                                                                      *
      SUBROUTINE PRWF(J)
!                                                                      *
!   Makes a (debug) printout of wave functions. There are two modes:   *
!                                                                      *
!      (1) J>0  - Used as a debug option in SOLVE, wavefunctions for   *
!                 orbital J are printed                                *
!      (2) J=0  - A  printout  of the grid and all wave functions is   *
!                 made                                                 *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP
      USE def_C
      USE grid_C
      USE int_C
      USE orb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE draw_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NB2, NROWS, II, II1, II2, K, MFK, JGG
      REAL(DOUBLE) :: CBZ
!-----------------------------------------------
!
!
      CBZ = C/Z
!
      IF (J > 0) THEN
!
!   Mode (1)
!
         WRITE (99, 300) NP(J), NH(J)
         NB2 = MTP0/2
         IF (2*NB2 == MTP0) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2 + 1
         ENDIF
         DO II = 1, NROWS
            II1 = II
            II2 = II1 + NROWS
            IF (II2 <= MTP0) THEN
               WRITE (99, 301) R(II1), P(II1), Q(II1), R(II2), P(II2), Q(II2)
            ELSE IF (II1 <= MTP0) THEN
               WRITE (99, 301) R(II1), P(II1), Q(II1)
            ENDIF
         END DO
         CALL DRAW (P, 1.0D00, Q, CBZ, MTP0)
!
      ELSE
!
!   Mode (2)
!
         DO K = 1, NW
            WRITE (99, 300) NP(K), NH(K)
            MFK = MF(K)
            NB2 = MFK/2
            IF (2*NB2 == MFK) THEN
               NROWS = NB2
            ELSE
               NROWS = NB2 + 1
            ENDIF
            DO II = 1, NROWS
               II1 = II
               II2 = II1 + NROWS
               IF (II2 <= MFK) THEN
                  WRITE (99, 301) R(II1), PF(II1,K), QF(II1,K), R(II2), PF(II2,&
                     K), QF(II2,K)
               ELSE IF (II1 <= MFK) THEN
                  WRITE (99, 301) R(II1), PF(II1,K), QF(II1,K)
               ENDIF
            END DO
            CALL DRAW (PF(:NNNP,K), 1.0D00, QF(:NNNP,K), CBZ, MF(K))
!
         END DO
!
      ENDIF
!
      RETURN
!
  300 FORMAT('1',1I2,1A2,' orbital:'/,2(&
         ' --------- r --------- ------- P (r) -------',&
         ' ------- Q (r) -------'))
  301 FORMAT(1P,6(1X,1D21.14))
      RETURN
!
      END SUBROUTINE PRWF
