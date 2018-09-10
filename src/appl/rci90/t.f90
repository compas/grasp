!***********************************************************************
! 
      function relci_qed_F_Klarsfeld(n,kappa,Z)                   result(F)
!--------------------------------------------------------------------
! Estimates the function  F (Z*\alpha) by using a series expansion
! from S Klarsfeld and A Maquet, Physics Letters  43B (1973) 201,
! Eqs (1) and (2) and the table of Bethe logarithms. The 
! vacuum-polarization contribution in Eq (2) is omitted. 
! This procedure is adapted from RCI92 of GRASP92, written
! by Farid A Parpia, to the Fortran 95 standard.
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
      !
      integer, intent(in)               :: n, kappa
      real(double), intent(in)         :: Z
      real(double)                     :: F
      !
      real(double), dimension(36), parameter :: bethe = &
         (/ 2.9841285_dp,   2.8117699_dp,  -0.0300167_dp,   2.7676636_dp, &
           -0.0381902_dp,  -0.0052321_dp,   2.7498118_dp,  -0.0419549_dp, &
           -0.0067409_dp,  -0.0017337_dp,   2.7408237_dp,  -0.0440347_dp, &
           -0.0076008_dp,  -0.0022022_dp,  -0.0007721_dp,   2.7356642_dp, &
           -0.0453122_dp,  -0.0081472_dp,  -0.0025022_dp,  -0.0009628_dp, &
           -0.0004079_dp,   2.7324291_dp,  -0.0461552_dp,  -0.0085192_dp, &
           -0.0027091_dp,  -0.0010945_dp,  -0.0004997_dp,  -0.0002409_dp, &
            2.7302673_dp,  -0.0467413_dp,  -0.0087850_dp,  -0.0028591_dp, &
           -0.0011904_dp,  -0.0005665_dp,  -0.0002904_dp,  -0.0001539_dp /) 
      !
      real(double), parameter :: C401 = 11.0_dp/24.0_dp,                 &
                               C402 = 3.0_dp/8.0_dp, ovlfac = 4.0_dp/3.0_dp
      !
      integer       :: l, loc
      real(double) :: bethel, factor, term
      !
      ! Ensure that the principal quantum number is in range
      if (n < 1   .or.   n > 8) then
         print *, "Principal quantum number,",n,", should be in the range 1-8."
         stop     "relci_qed_F_Klarsfeld(): program stop A."
      end if
      !
      l = angular_momentum_l(kappa)
      if (l > n-1) then
         print *, "Kappa = ",kappa," is out of range for n = ",n,"."
         stop     "relci_qed_F_Klarsfeld(): program stop B."
      end if
      !
      ! Find the appropriate entry in the table
      loc    = (n*n-n)/2+l+1
      bethel = bethe(loc)
      !
      ! Determine the quantity in square brackets in eq.(1) of
      ! Klarsfeld and Maquet
      term = -bethel
      !
      if (kappa > 0) then
         term = term - c402 / (l*(l+l+one))
      else
         term = term + c402 / ((l+one)*(l+l+one))
         if (kappa == -1) then
            factor = log (Z/c)
            factor = - (factor + factor)
            term   = term + factor + c401
         end if
      end if
      !
      F = ovlfac * term
      !
   end function relci_qed_F_Klarsfeld
