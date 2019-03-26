!***********************************************************************
!                                                                      *
    REAL(KIND(0.0D0))  FUNCTION VPINTF (IA,IB)
!                                                                      *
!   Computes nuclear vacuum polarization integrals.                    *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
     USE vast_kind_param, ONLY: DOUBLE
     USE parameter_def,   ONLY: KEYORB
     USE debug_C
     USE grid_C
     USE ncdist_C
     USE orb_C
     USE tatb_C, ONLY: mtp, ta, tb
     USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     USE quad_I
     IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IA, IB
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: k
!------------------------------------------------
!
      MTP = MIN (MF(IA),MF(IB))
      !TA(1) = 0.0D 00  !! Why would this be wrong?!!
      TA(1)  = 0
      DO 1 K = 2,MTP
         TA(K) = (PF(K,IA)*PF(K,IB)+QF(K,IA)*QF(K,IB))                  &
                 *ZDIST(K)
    1 CONTINUE
      CALL QUAD (VPINTF)
!
      IF (LDBPR(9)) WRITE (99,300) NP(IA),NH(IA),NP(IB),NH(IB),VPINTF
!
      RETURN
!
  300 FORMAT (/'VPINTF: V (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
!
      END FUNCTION VPINTF
