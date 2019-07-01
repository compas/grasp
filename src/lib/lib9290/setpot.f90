!***********************************************************************
!                                                                      *
      SUBROUTINE SETPOT(J, JP)
!                                                                      *
!   This  subroutine  sets  up the  arrays TF and TG for use by  the   *
!   subprograms IN, OUT, and SBSTEP.                                   *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      J:  (Input) Index of orbital                                    *
!      JP: (Output) Join point: point where TG changes sign            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:38   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE IOUNIT_C
      USE DEF_C
      USE GRID_C
      USE INT_C,           ONLY: TF, TG
      USE ORB_C
      USE POTE_C,          ONLY: YP
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: J
      INTEGER, INTENT(OUT) :: JP
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: DMHB2C, ENERGY, ENEFAC, RPI, YPRPOR
      LOGICAL :: NOTSET
!-----------------------------------------------
!
!
      NOTSET = .TRUE.
!
!   Define constants
!
      DMHB2C = -H/(2.0D00*C)
      ENERGY = E(J)
      ENEFAC = 2.0D00*C*C - ENERGY
!
!   Set up arrays TF and TG
!
!   Since TF(1) and TG(1) are never used, set them
!   to some arbitrary value
!
      JP = 0
      TF(1) = 0.0D00
      TG(1) = 0.0D00
      DO I = 2, N
         RPI = RP(I)
         YPRPOR = YP(I)*RPOR(I)
         TF(I) = DMHB2C*(ENEFAC*RPI + YPRPOR)
         TG(I) = DMHB2C*(ENERGY*RPI - YPRPOR)
         IF (.NOT.NOTSET) CYCLE
         IF (ABS(SIGN(1.0D00,TG(I))+SIGN(1.0D00,TG(1))) >= ACCY) CYCLE
         JP = I
         NOTSET = .FALSE.
      END DO
!
!   Trap for inappropriate grid
!
      IF (JP == 0) THEN
         WRITE (ISTDE, *) 'SETPOT: Join set to NNNP/2 = ', nnnp/2
         JP = NNNP/2
      ENDIF
!
      RETURN
      END SUBROUTINE SETPOT
