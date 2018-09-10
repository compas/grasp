!***********************************************************************
!                                                                      *
      SUBROUTINE RKINTC(IA, IB, IC, ID, K, TEGRAL) 
!                                                                      *
!                         k                                            *
!   This routine returns R (abcd) integrals.                           *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE cteilsrk_C
      USE orb_C
      USE kkstart_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: IA, IB, IC, ID 
      INTEGER, INTENT(IN) :: K 
      REAL(DOUBLE), INTENT(OUT) :: TEGRAL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KEY, ISWAP, INDEX, JL, JU, JM, LOC 
      LOGICAL :: FOUND, FIRST 
!-----------------------------------------------
!
      KEY = NW + 1 
!
!   Ensure that the indices are in `canonical' order
!   Compute the composite (packed) index
!
 
 
      IF (IA > IC) THEN 
         ISWAP = IC 
         IC = IA 
         IA = ISWAP 
      ENDIF 
      IF (IB > ID) THEN 
         ISWAP = ID 
         ID = IB 
         IB = ISWAP 
      ENDIF 
      IF (IA > IB) THEN 
         ISWAP = IB 
         IB = IA 
         IA = ISWAP 
         ISWAP = ID 
         ID = IC 
         IC = ISWAP 
      ENDIF 
 
      INDEX = ((IA*KEY + IB)*KEY + IC)*KEY + ID 
!
      JL = KSTART(K) 
      JU = KSTART(K+1) - 1 
 
      IF (INDEX<INDTEIRK(JL) .OR. INDEX>INDTEIRK(JU)) THEN 
         WRITE (*, *) 'Something wrong in rkintc' 
         STOP  
      ENDIF 
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 CONTINUE 
      IF (JU - JL > 1) THEN 
         JM = (JU + JL)/2 
         IF (INDTEIRK(JM) > INDEX) THEN 
            JU = JM 
         ELSE 
            JL = JM 
         ENDIF 
         GO TO 1 
      ENDIF 
!
!   The range is bracketed to the extent possible
!
      IF (INDEX == INDTEIRK(JU)) THEN 
         LOC = JU 
      ELSE IF (INDEX == INDTEIRK(JL)) THEN 
         LOC = JL 
      ELSE 
         WRITE (*, *) K, IA, IB, IC, ID, INDEX 
         WRITE (*, *) 'Rkintc Integral not found' 
         STOP  
      ENDIF 
!
!   Return the value of the integral
!   from storage
 
      TEGRAL = VALTEIRK(LOC) 
!
      RETURN  
      END SUBROUTINE RKINTC 
