!***********************************************************************
!                                                                      *
      SUBROUTINE ENGOUT(EAV, E, JTOT, IPAR, ILEV, NN, MODE)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:24   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEF_C, ONLY: AUCM, AUEV
      USE JLABL_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NN
      INTEGER, INTENT(IN) :: MODE
      REAL(DOUBLE), INTENT(IN) :: EAV
      INTEGER, INTENT(IN) :: JTOT(NN)
      INTEGER, INTENT(IN) :: IPAR(NN)
      INTEGER, INTENT(IN) :: ILEV(NN)
      REAL(DOUBLE), INTENT(IN) :: E(NN)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: J, I, IP
      REAL(DOUBLE) :: EAU, ECM, EEV
!-----------------------------------------------
!
!   Always print the eigenenergies
!
      WRITE (24, 300)
      WRITE (24, 301)
      DO J = 1, NN
         I = ILEV(J)
         EAU = E(J) + EAV
         ECM = EAU*AUCM
         EEV = EAU*AUEV
         IP = (IPAR(J)+3)/2
         WRITE (24, 302) I, JLBR(JTOT(J)), JLBP(IP), EAU, ECM, EEV
!          WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
      END DO
!
      IF (NN > 1) THEN
!
!   Energy separations
!
         IF (MODE==1 .OR. MODE==3) THEN
            WRITE (24, 303)
            WRITE (24, 301)
            DO J = 2, NN
               I = ILEV(J)
               EAU = E(J) - E(J-1)
               ECM = EAU*AUCM
               EEV = EAU*AUEV
               IP = (IPAR(J)+3)/2
               WRITE (24, 302) I, JLBR(JTOT(J)), JLBP(IP), EAU, ECM, EEV
!               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
            END DO
         ENDIF
!
!   Energies relative to level 1
!
         IF (MODE==2 .OR. MODE==3) THEN
            WRITE (24, 304)
            WRITE (24, 301)
            DO J = 2, NN
               I = ILEV(J)
               EAU = E(J) - E(1)
               ECM = EAU*AUCM
               EEV = EAU*AUEV
               IP = (IPAR(J)+3)/2
               WRITE (24, 302) I, JLBR(JTOT(J)), JLBP(IP), EAU, ECM, EEV
!               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
            END DO
         ENDIF
!
      ENDIF
!
      RETURN
!
  300 FORMAT(/,'Eigenenergies:')
  301 FORMAT(/,'Level  J Parity',7X,'Hartrees',14X,'Kaysers',16X,'eV'/)
  302 FORMAT(1I3,2X,2A4,1P,3D22.14)
  303 FORMAT(/,'Energy of each level relative to immediately lower',' level:')
  304 FORMAT(/,'Energy of each level relative to lowest level:')
      RETURN
!
      END SUBROUTINE ENGOUT
