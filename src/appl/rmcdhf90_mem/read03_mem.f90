!*******************************************************************
!                                                                  *
      SUBROUTINE read03_mem(JBLOCK)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE hmat_C,     ONLY: IROW, IENDC
      USE rmcdhf_mem_C
      USE mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: JBLOCK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------

      DO I = IMIN1_30(JBLOCK)+1+MYID, IMAX1_30(JBLOCK), NPROCS
         IENDC(I-IMIN1_30(JBLOCK)) = IENDC_30(I)
      END DO
      DO I = IMIN2_30(JBLOCK)+1, IMAX2_30(JBLOCK)
         IROW(I-IMIN2_30(JBLOCK)) = IROW_30(I)
      END DO

      RETURN
      END SUBROUTINE read03_mem
