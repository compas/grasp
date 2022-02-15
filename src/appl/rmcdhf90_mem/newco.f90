!***********************************************************************
!                                                                      *
      SUBROUTINE NEWCO(SUM)
!                                                                      *
!   This routine computes the level weights, the generalized occupa-   *
!   tion numbers,  and average energy for  EOL calculations;  this     *
!   information and the eigenvectors are then printed out.             *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!               [RSCF92]: CSFWGT, DSUBRS.                              *
!                                         Last revision: 24 Dec 1992   *
!   Block version by Xinghong He          Last revision: 05 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:53:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE DEBUG_C
      USE def_C
      USE eigv_C
      USE orb_C
      USE scf_C
      USE syma_C,          ONLY: JPGG
      USE iounit_C
      USE hblock_C
      USE pos_C
      USE peav_C
      USE blkidx_C
      USE mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dsubrs_I
      USE iq_I
      USE csfwgt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(OUT) :: SUM
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, NB, I, NOFF, JALL, JBLOCK
      REAL(DOUBLE) :: WEITJ, EE
      LOGICAL :: EOL
!-----------------------------------------------
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!
!   Compute weighting factors
!
      SUM = 0.D0
      DO J = 1, NCMIN
         WEITJ = WEIGHT(J)
         IF (WEITJ == (-2.D0)) THEN
            WT(J) = 1.D0
         ELSE IF (WEITJ == (-1.D0)) THEN
!GGGG
!GG            WT(J) = IATJPO(J)
            jblock = idxblk(J)     ! Block number of this state
            WT(J) = ABS(JPGG(jblock))
         ELSE
            WT(J) = WEIGHT(J)
         ENDIF
         SUM = SUM + WT(J)
      END DO

      WT(:NCMIN) = WT(:NCMIN)/SUM
!
!   Compute generalised occupation numbers
!   <----- Distributed ----->
!
      EOL = .TRUE.

      DO J = 1, NW
         SUM = 0.D0
         DO NB = 1, NBLOCK
            DO I = MYID + 1, NCFBLK(NB), NPROCS
               SUM = SUM + DSUBRS(EOL,I,I,NB)*IQ(J,I + NCFPAST(NB))
            END DO
         END DO
         UCF(J) = SUM
      END DO
!
!   Write out level energies and weights
!
      WRITE (*, 300)
      SUM = 0.D0
      NOFF = 0
      DO JALL = 1, NCMIN
         NB = IDXBLK(JALL)                       ! Block number of this state
         EE = EAVBLK(NB) + EVAL(JALL)
         WRITE (*, 301) ICCMIN(JALL), EE, WT(JALL)
         IF (LDBPG(5)) THEN
            WRITE (99, *) JALL, NB, NCFBLK(NB), NEVECPAST(NB)
            WRITE (99, 302)
            WRITE (99, 303) (EVEC(I + NOFF),I=1,NCFBLK(NB))
            NOFF = NOFF + NCFBLK(NB)
         ENDIF
         SUM = SUM + WT(JALL)*EE
      END DO

      CALL CSFWGT (.TRUE.)
!
!   Write out average energy
!
      IF (NCMIN > 1) WRITE (*, 304) SUM
!
!   Write out generalized occupation numbers
!
      WRITE (*, 305)
      WRITE (*, 303) (UCF(I),I=1,NW)

  300 FORMAT(/,'Optimise on the following level(s):'/)
  301 FORMAT('Level ',1I2,4X,'Energy = ',1P,1D19.12,4X,'Weight = ',1D12.5)
  302 FORMAT(/,'Configuration mixing coefficients:')
  303 FORMAT(1X,1P,6D12.4)
  304 FORMAT(/,'Weighted average energy of these levels = ',1P,D18.10)
  305 FORMAT(/,'Generalised occupation numbers:'/)

      RETURN
      END SUBROUTINE NEWCO
