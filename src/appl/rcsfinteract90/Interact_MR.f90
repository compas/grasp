!***********************************************************************
!                                                                      *
      SUBROUTINE Interact_MR(ICOLBREI,I_Count)
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                       k                              *
!                   T  (ab)            V  (abcd)                       *
!                    rs                 rs                             *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   c and d are orbital sequence numbers.  r and s are configuration   *
!   state function indices.                                            *
!                                                                      *
!   Written by  G. Gaigalas                      NIST, December 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE BLK_C,            only:  NBLOCK,NCFBLK
      USE rang_Int_C
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: ICOLBREI
      INTEGER, INTENT(INOUT) :: I_Count
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JA, JB, int_CSF
!-----------------------------------------------
!
      DO JA = 1, NCFBLK(NBLOCK)
         JB = NCFBLK(NBLOCK) + 1
         CALL Interact_CSF(JA,JB,ICOLBREI,int_CSF)
         if(int_CSF .NE. 0) THEN
            I_count = I_count + 1
            WRITE(22,'(A)') TRIM(C_shell(JB))
            WRITE(22,'(A)') TRIM(C_quant(JB))
            WRITE(22,'(A)') TRIM(C_coupl(JB))
            RETURN
         END IF
      END DO
      RETURN
      END SUBROUTINE Interact_MR
