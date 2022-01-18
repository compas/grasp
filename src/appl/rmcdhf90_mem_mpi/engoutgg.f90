!***********************************************************************
!                                                                      *
      SUBROUTINE ENGOUTGG(E,ILEV,NN,MODE)
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
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE def_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE blkidx_C
      USE peav_C
      USE syma_C,          ONLY: JPGG
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NN
      INTEGER, INTENT(IN) :: MODE
      INTEGER, INTENT(IN) :: ILEV(NN)
      REAL(DOUBLE), INTENT(IN) :: E(NN)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, JBLOCK, I, IP, JTOT
      REAL(DOUBLE) :: EAV, EAU, ECM, EEV
!-----------------------------------------------
!
!
!   Always print the eigenenergies
!
      WRITE (24, 300)
      WRITE (24, 301)
      DO J = 1, NN
         JBLOCK = IDXBLK(J)
         EAV = EAVBLK(JBLOCK)
         I = ILEV(J)
         EAU = E(J) + EAV
         ECM = EAU*AUCM
         EEV = EAU*AUEV
!GG         IP = (IPAR(J)+3)/2
         JTOT = IABS(JPGG(jblock))
         IF(JPGG(jblock) >= 0) THEN
            ip = 2
         ELSE
            ip = 1
         END IF
         WRITE (24, 302) I, LABJ(JTOT), LABP(IP), EAU, ECM, EEV
!GG         WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
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
!GG               IP = (IPAR(J)+3)/2
               jblock = idxblk(j)
               JTOT = IABS(JPGG(jblock))
               IF(JPGG(jblock) >= 0) THEN
                  ip = 2
               ELSE
                  ip = 1
               END IF
               WRITE (24, 302) I, LABJ(JTOT), LABP(IP), EAU, ECM, EEV
!GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
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
!GG               IP = (IPAR(J)+3)/2
               jblock = idxblk(j)
               JTOT = IABS(JPGG(jblock))
               IF(JPGG(jblock) >= 0) THEN
                  ip = 2
               ELSE
                  ip = 1
               END IF
               WRITE (24, 302) I, LABJ(JTOT), LABP(IP), EAU, ECM, EEV
!GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
            END DO
         ENDIF
!
      ENDIF
!
      RETURN
!
  300 FORMAT(/,'Eigenenergies:')
  301 FORMAT(/,'Level  J Parity',7X,'Hartrees',14X,'Kaysers',16X,'eV'/)
  302 FORMAT(1I3,4X,2A4,1P,3D21.12)
  303 FORMAT(/,'Energy of each level relative to immediately lower',' level:')
  304 FORMAT(/,'Energy of each level relative to lowest level:')
      RETURN
!
      END SUBROUTINE ENGOUTGG
