!***********************************************************************
!                                                                      *
      SUBROUTINE VACPOL
!                                                                      *
!   This  routine controls the setting up of the vacuum polarization   *
!   potential for the given nuclear charge distribution at each grid   *
!   point using the analytic functions defined by  L Wayne Fullerton   *
!   and G A Rinker Jr in Phys Rev A 13 (1976) 1283. The potential is   *
!   accumulated  in  array TB(I), I = 1, ..., N  which  is in COMMON   *
!   block /TATB/ .                                                     *
!                                                                      *
!   Call(s) to: [RCI92]: VAC2, VAC4.                                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE grid_C
      USE npar_C
      USE ncdist_C
      USE tatb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE vac2_I
      USE vac4_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
!-----------------------------------------------
!
!   Redefine  ZDIST  to be  rho*r*r'
!
      ZDIST(:MTP) = ZDIST(:MTP)*R(:MTP)*RP(:MTP)
!
!   Second-order vacuum polarisation potential; returned in
!   array TB
!
      CALL VAC2
!
!   Fourth-order vacuum polarization potential; returned in
!   array TA
!
      CALL VAC4
!
!   If option 7 is set, use user-defined vacuum polarization
!   potential
!
      RETURN
      END SUBROUTINE VACPOL
