!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION DSUBRS (EOL, I, J, JBLOCK) 
!                                                                      *
!   The coefficients d   for  I = r, J = s  are calculated here.       *
!                     rs                                               *
!                                                                      *
!                         NCMIN                                        *
!                          Sum     w     c          c                  *
!                         t = 1     t     r Gamma    s Gamma           *
!                                                t          t          *
!    EOL:          d   = -------------------------------------         *
!                   rs             NCMIN                               *
!                                   Sum     w                          *
!                                  t = 1     t                         *
!                                                                      *
!                                                                      *
!                                           /  NCF                     *
!    EAL:             d   =  delta     w   /   Sum   w                 *
!                      rs         rs    r /   t = 1   t                *
!                                                                      *
!                                                                      *
!   Written by Farid A Parpia               Last update: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!   This routine assumes (and does not check) that I and J belong
!   the same block and the job is done within the block specified
!   by parameter jblock
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE def_C
      USE eigv_C
      USE hblock_C
      USE pos_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER , INTENT(IN) :: J 
      INTEGER  :: JBLOCK 
      LOGICAL , INTENT(IN) :: EOL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KCMIN, IK, JK 
!-----------------------------------------------
 
      IF (EOL) THEN 
         DSUBRS = 0.D0 
         DO K = 1, NEVBLK(JBLOCK) 
            KCMIN = K + NCMINPAST(JBLOCK) 
            IK = I + NEVECPAST(JBLOCK) + (K - 1)*NCFBLK(JBLOCK) 
            JK = J + NEVECPAST(JBLOCK) + (K - 1)*NCFBLK(JBLOCK) 
            DSUBRS = DSUBRS + EVEC(IK)*EVEC(JK)*WT(KCMIN) 
         END DO 
      ELSE IF (I == J) THEN 
         DSUBRS = WT(I) 
      ELSE 
         DSUBRS = 0.D0 
      ENDIF 
 
      RETURN  
      END FUNCTION DSUBRS 
