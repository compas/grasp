!***********************************************************************
!                                                                      *
      SUBROUTINE YPOT(J) 
!                                                                      *
!   This subroutine  tabulates the  potential function Y(r) (Eq (14)   *
!   in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C   *
!   Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The   *
!   function is tabulated in the COMMON array YP.                      *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK                                     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: KEYORB
      USE debug_C
      USE grid_C
      USE npot_C
      USE orb_C
      USE pote_C
      USE scf_C
      USE tatb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE yzk_I 
      USE draw_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: J 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INDEX, LABEL, K, IOY1, IOY2, NB3, NROWS, II, II1, II2, II3 
      REAL(DOUBLE) :: COEFF 
!-----------------------------------------------
!   Debug printout: composition
 
      IF (LDBPR(29) .OR. LDBPR(30)) WRITE (99, 300) NP(J), NH(J) 
!
!   Initialize array YP with the nuclear potential piece
!
!   Since YA() below  contains contributions from THIS node only,
!   the initialization should be in consistence with that.
 
      YP(:N) = ZZ(:N)
 
      DO INDEX = 1, NYCOF 
 
         ! Decode information in label
         LABEL = NYA(INDEX) 
         K = MOD(LABEL,KEY) 
         LABEL = LABEL/KEY 
         IOY1 = MOD(LABEL,KEY) 
         IOY2 = LABEL/KEY 
         COEFF = YA(INDEX) 
 
         IF (LDBPR(29)) WRITE (99, 301) K, COEFF, NP(IOY1), NH(IOY1), NP(IOY2)&
            , NH(IOY2) 
 
         CALL YZK (K, IOY1, IOY2)                ! Accumulate contributions 
         YP(:N) = YP(:N) - COEFF*TB(:N) 
      END DO 
!
!   Debug printout
!
      IF (LDBPR(30)) THEN 
         WRITE (99, 302) 
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
               WRITE (99, 303) R(II1), YP(II1), R(II2), YP(II2), R(II3), YP(II3&
                  ) 
            ELSE IF (II2 <= N) THEN 
               WRITE (99, 303) R(II1), YP(II1), R(II2), YP(II2) 
            ELSE 
               WRITE (99, 303) R(II1), YP(II1) 
            ENDIF 
         END DO 
         CALL DRAW (YP, 1.0D00, YP, 0.0D00, N) 
      ENDIF 
!
      RETURN  
!
  300 FORMAT(/,/,' Direct potential for ',1I2,1A2,' orbital :'/,/) 
  301 FORMAT(/,25X,'(',1I2,')'/,1X,1P,D21.14,'* Y    (',1I2,1A2,',',1I2,1A2,')'&
         ) 
  302 FORMAT(/,/,3(' --------- r --------- ------- Y (r) -------')) 
  303 FORMAT(1P,6(1X,1D21.14)) 
      RETURN  
!
      END SUBROUTINE YPOT 
