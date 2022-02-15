!*******************************************************************
!                                                                  *
      SUBROUTINE read1_mem(NFILE, NCONTR_tot2, LAB, NCONTR)
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
      INTEGER        :: NFILE, LAB, NCONTR
      INTEGER(LONG)  :: NCONTR_tot2
!-----------------------------------------------

      if(NFILE == 31) then
         LAB    = LAB_31(NCONTR_tot2)
         NCONTR = NCONTR_31(NCONTR_tot2)
      else if(NFILE == 32) then
         LAB    = LAB_32(NCONTR_tot2)
         NCONTR = NCONTR_32(NCONTR_tot2)
      else if(NFILE == 33) then
         LAB    = LAB_33(NCONTR_tot2)
         NCONTR = NCONTR_33(NCONTR_tot2)
      else if(NFILE == 34) then
         LAB    = LAB_34(NCONTR_tot2)
         NCONTR = NCONTR_34(NCONTR_tot2)
      else if(NFILE == 35) then
         LAB    = LAB_35(NCONTR_tot2)
         NCONTR = NCONTR_35(NCONTR_tot2)
      else if(NFILE == 36) then
         LAB    = LAB_36(NCONTR_tot2)
         NCONTR = NCONTR_36(NCONTR_tot2)
      else if(NFILE == 37) then
         LAB    = LAB_37(NCONTR_tot2)
         NCONTR = NCONTR_37(NCONTR_tot2)
      else if(NFILE == 38) then
         LAB    = LAB_38(NCONTR_tot2)
         NCONTR = NCONTR_38(NCONTR_tot2)
      else if(NFILE == 39) then
         LAB    = LAB_39(NCONTR_tot2)
         NCONTR = NCONTR_39(NCONTR_tot2)
      else if(NFILE == 40) then
         LAB    = LAB_40(NCONTR_tot2)
         NCONTR = NCONTR_40(NCONTR_tot2)
      else if(NFILE == 41) then
         LAB    = LAB_41(NCONTR_tot2)
         NCONTR = NCONTR_41(NCONTR_tot2)
      else if(NFILE == 42) then
         LAB    = LAB_42(NCONTR_tot2)
         NCONTR = NCONTR_42(NCONTR_tot2)
      else if(NFILE == 43) then
         LAB    = LAB_43(NCONTR_tot2)
         NCONTR = NCONTR_43(NCONTR_tot2)
      else if(NFILE == 44) then
         LAB    = LAB_44(NCONTR_tot2)
         NCONTR = NCONTR_44(NCONTR_tot2)
      else if(NFILE == 45) then
         LAB    = LAB_45(NCONTR_tot2)
         NCONTR = NCONTR_45(NCONTR_tot2)
      else if(NFILE == 46) then
         LAB    = LAB_46(NCONTR_tot2)
         NCONTR = NCONTR_46(NCONTR_tot2)
      else if(NFILE == 47) then
         LAB    = LAB_47(NCONTR_tot2)
         NCONTR = NCONTR_47(NCONTR_tot2)
      else if(NFILE == 48) then
         LAB    = LAB_48(NCONTR_tot2)
         NCONTR = NCONTR_48(NCONTR_tot2)
      else if(NFILE == 49) then
         LAB    = LAB_49(NCONTR_tot2)
         NCONTR = NCONTR_49(NCONTR_tot2)
      else if(NFILE == 50) then
         LAB    = LAB_50(NCONTR_tot2)
         NCONTR = NCONTR_50(NCONTR_tot2)
      else
         print*, "Error in read1_memi  NFILE=",NFILE
         stop
      end if

      RETURN
      END SUBROUTINE read1_mem
