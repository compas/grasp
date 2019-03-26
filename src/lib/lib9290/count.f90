!***********************************************************************
!                                                                      *
      SUBROUTINE COUNT(FR, MTPFR, NNCFF, SGN)
!                                                                      *
!   This subroutine counts the nodes in the radial function FR using   *
!   the criteria  given by C Froese Fischer, Comp Phys Rep, 3 (1986)   *
!   314-315 . The  function FR is assumed defined on the first MTPFR   *
!   points of the radial grid. The sign of the function at the first   *
!   oscillation is also determined.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:58   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------                       *
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP
      USE COUN_C
      USE DEF_C, ONLY: ACCY
      USE GRID_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: MTPFR
      INTEGER, INTENT(OUT) :: NNCFF
      REAL(DOUBLE), INTENT(OUT) :: SGN
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: FR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(NNNP) :: LCEXT
      INTEGER :: NEXT, I, LOC, NSTPS
      REAL(DOUBLE) :: EXT, EMX, ABFRI, TEST, THRESE, ABLCL
!-----------------------------------------------
!
!
!   (1) Find all extrema in FR
!   (2) Find the maximum amplitudes of FR
!
      NEXT = 1
      EXT = 0.0D00
      LCEXT(1) = 1
      EMX = 0.0D00
      DO I = 2, MTPFR
         ABFRI = ABS(FR(I))
         TEST = ABS(SIGN(1.0D00,FR(I))+SIGN(1.0D00,FR(I-1)))
         IF (TEST <= ACCY) THEN
            NEXT = NEXT + 1
            LCEXT(NEXT) = 0
            EXT = 0.0D00
         ENDIF
         IF (ABFRI > EXT) THEN
            EXT = ABFRI
            LCEXT(NEXT) = I
         ENDIF
         IF (ABFRI <= EMX) CYCLE
         EMX = ABFRI
      END DO
!
!   Eliminate oscillations with amplitude less than THRESH times
!   the maximum
!
      LOC = 0
      THRESE = THRESH*EMX
    4 CONTINUE
      LOC = LOC + 1
      IF (LOC <= NEXT) THEN
         IF (LCEXT(LOC) == 0) THEN
            ABLCL = 0.0D00
         ELSE
            ABLCL = ABS(FR(LCEXT(LOC)))
         ENDIF
         IF (ABLCL < THRESE) THEN
            NEXT = NEXT - 1
            NSTPS = NEXT - LOC
            LCEXT(LOC:NSTPS+LOC) = LCEXT(LOC+1:NSTPS+1+LOC)
            LOC = LOC - 1
         ENDIF
         GO TO 4
      ENDIF
!
!   Count changes of sign using the remaining oscillations
!
      NNCFF = 0
      DO I = 2, NEXT
         TEST = ABS(SIGN(1.0D00,FR(LCEXT(I)))+SIGN(1.0D00,FR(LCEXT(I-1))))
         IF (TEST > ACCY) CYCLE
         NNCFF = NNCFF + 1
      END DO
!
!   Determine the position of the first oscillation, and the
!   sign of the function at this location
!
      SGN = SIGN(1.0D00,FR(LCEXT(1)))
!
      RETURN
      END SUBROUTINE COUNT
