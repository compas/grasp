!***********************************************************************
!                                                                      *
      SUBROUTINE ENGOUT1(EAV, E, JTOT, IPAR, ILEV, NN, MODE, K) 
!                                                                      *
!   This  subroutine prints  energy levels, splittings, and energies   *
!   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
!   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
!   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
!   and separations are printed. If MODE is 2, the eigenenergies and   *
!   energies relative to level 1 are printed. If MODE is 3, the eig-   *
!   enenergies,  separations,  and energies relative to level  1 are   *
!   printed.                                                           *
!                                          Last updated: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE def_C, ONLY: AUCM, AUEV, CCMS, FASI, FBSI
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NN 
      INTEGER  :: MODE 
      INTEGER , INTENT(IN) :: K 
      REAL(DOUBLE) , INTENT(IN) :: EAV 
      INTEGER , INTENT(IN) :: JTOT(NN) 
      INTEGER , INTENT(IN) :: IPAR(NN) 
      INTEGER , INTENT(IN) :: ILEV(NN) 
      REAL(DOUBLE) , INTENT(IN) :: E(NN) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, IP 
      REAL(DOUBLE) :: EAU, ECM, EEV 
!-----------------------------------------------
!
!   Always print the eigenenergies
!
      IF (K == 1) WRITE (24, 299) 
      IF (K == 2) WRITE (24, 300) 
      WRITE (24, 301) 
      DO J = 1, NN 
         I = ILEV(J) 
         EAU = E(J) + EAV 
         ECM = EAU*AUCM 
         EEV = EAU*AUEV 
         IP = (IPAR(J)+3)/2 
         WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM 
      END DO 
!
      RETURN  
!
  299 FORMAT('Eigenenergies for the initial state list') 
  300 FORMAT('Eigenenergies for the final state list') 
  301 FORMAT('Level  J Parity',10X,'Hartrees',18X,'Kaysers') 
  302 FORMAT(1I3,4X,2A4,1P,2D25.15) 
  303 FORMAT('Energy of each level relative to immediately lower',' level:') 
  304 FORMAT('Energy of each level relative to lowest level:') 
      RETURN  
!
      END SUBROUTINE ENGOUT1 
