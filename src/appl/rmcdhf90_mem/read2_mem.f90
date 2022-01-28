!*******************************************************************
!                                                                  *
      SUBROUTINE read2_mem(NFILE,NCONTR_tot,NCONTR,ICLMNDUM)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: LONG
      USE rmcdhf_mem_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER(LONG) :: NCONTR_tot
      INTEGER       :: NFILE, NCONTR, ICLMNDUM
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------

      DO I = 1, NCONTR
      !*** Read column index, sparse matrix index, and
      !*** coefficient for all contributions from this integral
         if(NFILE == 31) then
            ICLMNDUM = ICLMN_31(NCONTR_tot+I)
            INDX(I) = INDX_31(NCONTR_tot+I)
            COEFF(I) = COEFF_31(NCONTR_tot+I)
         else if(NFILE == 32) then
            ICLMNDUM = ICLMN_32(NCONTR_tot+I)
            INDX(I) = INDX_32(NCONTR_tot+I)
            COEFF(I) = COEFF_32(NCONTR_tot+I)
         else if(NFILE == 33) then
            ICLMNDUM = ICLMN_33(NCONTR_tot+I)
            INDX(I) = INDX_33(NCONTR_tot+I)
            COEFF(I) = COEFF_33(NCONTR_tot+I)
         else if(NFILE == 34) then
            ICLMNDUM = ICLMN_34(NCONTR_tot+I)
            INDX(I) = INDX_34(NCONTR_tot+I)
            COEFF(I) = COEFF_34(NCONTR_tot+I)
         else if(NFILE == 35) then
            ICLMNDUM = ICLMN_35(NCONTR_tot+I)
            INDX(I) = INDX_35(NCONTR_tot+I)
            COEFF(I) = COEFF_35(NCONTR_tot+I)
         else if(NFILE == 36) then
            ICLMNDUM = ICLMN_36(NCONTR_tot+I)
            INDX(I) = INDX_36(NCONTR_tot+I)
            COEFF(I) = COEFF_36(NCONTR_tot+I)
         else if(NFILE == 37) then
            ICLMNDUM = ICLMN_37(NCONTR_tot+I)
            INDX(I) = INDX_37(NCONTR_tot+I)
            COEFF(I) = COEFF_37(NCONTR_tot+I)
         else if(NFILE == 38) then
            ICLMNDUM = ICLMN_38(NCONTR_tot+I)
            INDX(I) = INDX_38(NCONTR_tot+I)
            COEFF(I) = COEFF_38(NCONTR_tot+I)
         else if(NFILE == 39) then
            ICLMNDUM = ICLMN_39(NCONTR_tot+I)
            INDX(I) = INDX_39(NCONTR_tot+I)
            COEFF(I) = COEFF_39(NCONTR_tot+I)
         else if(NFILE == 40) then
            ICLMNDUM = ICLMN_40(NCONTR_tot+I)
            INDX(I) = INDX_40(NCONTR_tot+I)
            COEFF(I) = COEFF_40(NCONTR_tot+I)
         else if(NFILE == 41) then
            ICLMNDUM = ICLMN_41(NCONTR_tot+I)
            INDX(I) = INDX_41(NCONTR_tot+I)
            COEFF(I) = COEFF_41(NCONTR_tot+I)
         else if(NFILE == 42) then
            ICLMNDUM = ICLMN_42(NCONTR_tot+I)
            INDX(I) = INDX_42(NCONTR_tot+I)
            COEFF(I) = COEFF_42(NCONTR_tot+I)
         else if(NFILE == 43) then
            ICLMNDUM = ICLMN_43(NCONTR_tot+I)
            INDX(I) = INDX_43(NCONTR_tot+I)
            COEFF(I) = COEFF_43(NCONTR_tot+I)
         else if(NFILE == 44) then
            ICLMNDUM = ICLMN_44(NCONTR_tot+I)
            INDX(I) = INDX_44(NCONTR_tot+I)
            COEFF(I) = COEFF_44(NCONTR_tot+I)
         else if(NFILE == 45) then
            ICLMNDUM = ICLMN_45(NCONTR_tot+I)
            INDX(I) = INDX_45(NCONTR_tot+I)
            COEFF(I) = COEFF_45(NCONTR_tot+I)
         else if(NFILE == 46) then
            ICLMNDUM = ICLMN_46(NCONTR_tot+I)
            INDX(I) = INDX_46(NCONTR_tot+I)
            COEFF(I) = COEFF_46(NCONTR_tot+I)
         else if(NFILE == 47) then
            ICLMNDUM = ICLMN_47(NCONTR_tot+I)
            INDX(I) = INDX_47(NCONTR_tot+I)
            COEFF(I) = COEFF_47(NCONTR_tot+I)
         else if(NFILE == 48) then
            ICLMNDUM = ICLMN_48(NCONTR_tot+I)
            INDX(I) = INDX_48(NCONTR_tot+I)
            COEFF(I) = COEFF_48(NCONTR_tot+I)
         else if(NFILE == 49) then
            ICLMNDUM = ICLMN_49(NCONTR_tot+I)
            INDX(I) = INDX_49(NCONTR_tot+I)
            COEFF(I) = COEFF_49(NCONTR_tot+I)
         else if(NFILE == 50) then
            ICLMNDUM = ICLMN_50(NCONTR_tot+I)
            INDX(I) = INDX_50(NCONTR_tot+I)
            COEFF(I) = COEFF_50(NCONTR_tot+I)
         else
            print*, "Error in read2_mem"
            stop
         end if
      END DO

      RETURN
      END SUBROUTINE read2_mem
