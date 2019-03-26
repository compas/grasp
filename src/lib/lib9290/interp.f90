!***********************************************************************
!                                                                      *
      SUBROUTINE INTERP(XARR, YARR, NARR, XVAL, YVAL, ACCY)
!                                                                      *
!   This routine returns  YVAL  given a value  XVAL by interpolating   *
!   using a pair of arrays XARR(1:NARR), YARR(1:NARR), that tabulate   *
!   a  function.  ACCY  is the  desired  accuracy of the estimate: a   *
!   warning message is issued if this is not achieved.  A warning is   *
!   also issued when the routine is extrapolating.  Aitken's algori-   *
!   thm is used. See, for  instance, F B Hildebrand, Introduction to   *
!   Numerical Analysis,  2nd ed., McGraw-Hill, New York, NY, 1974.     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:32   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NARR
      REAL(DOUBLE), INTENT(IN) :: XVAL
      REAL(DOUBLE), INTENT(OUT) :: YVAL
      REAL(DOUBLE), INTENT(IN) :: ACCY
      REAL(DOUBLE), DIMENSION(NARR), INTENT(IN) :: XARR
      REAL(DOUBLE), DIMENSION(NARR), INTENT(IN) :: YARR
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: MXORD = 11
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NRSTLO, NRSTHI, K, LLO, LHI, LLR, IROW, LOCNXT, ILIROK, ILDIAG&
         , ILOTHR, IBEST
      REAL(DOUBLE), DIMENSION(MXORD) :: DX, X, EST
      REAL(DOUBLE), DIMENSION((MXORD*(MXORD + 1))/2) :: POLY
      REAL(DOUBLE) :: DIFF, DIFFT, DEBEB, DEBE
      LOGICAL :: SET
      LOGICAL, DIMENSION(2*MXORD + 2) :: USED
!-----------------------------------------------
!
!
!   MXORD is the maximum order of the interpolation
!
!
!   This function dispenses with the need for a two-dimensional
!   array for the Aitken lower triangle matrix
!
!
!   Determine the nearest two XARR entries bounding XVAL
!
      IF (XVAL < XARR(1)) THEN
         NRSTLO = 1
         NRSTHI = 1
         WRITE (*, 300)
      ELSE IF (XVAL > XARR(NARR)) THEN
         NRSTLO = NARR
         NRSTHI = NARR
         WRITE (*, 300)
      ELSE
         K = 0
    1    CONTINUE
         K = K + 1
         IF (XARR(K) < XVAL) THEN
            NRSTLO = K
            GO TO 1
         ELSE
            NRSTHI = K
         ENDIF
      ENDIF
!
!   Clear relevant piece of use-indicator array
!
      LLO = MAX(NRSTLO - MXORD,1)
      LHI = MIN(NRSTHI + MXORD,NARR)
      LLR = LLO - 1
      USED(LLO-LLR:LHI-LLR) = .FALSE.
!
!   Determine next-nearest XARR entry
!
      DO IROW = 1, MXORD
         LLO = MAX(NRSTLO - IROW + 1,1)
         LHI = MIN(NRSTHI + IROW - 1,NARR)
         SET = .FALSE.
         DO K = LLO, LHI
            IF (USED(K-LLR)) CYCLE
            IF (.NOT.SET) THEN
               DIFF = XARR(K) - XVAL
               LOCNXT = K
               SET = .TRUE.
            ELSE
               DIFFT = XARR(K) - XVAL
               IF (ABS(DIFFT) < ABS(DIFF)) THEN
                  DIFF = DIFFT
                  LOCNXT = K
               ENDIF
            ENDIF
         END DO
         USED(LOCNXT-LLR) = .TRUE.
         X(IROW) = XARR(LOCNXT)
         DX(IROW) = DIFF
!
!   Fill table for this row
!
         DO K = 1, IROW
            ILIROK = ILOC(IROW,K)
            IF (K == 1) THEN
               POLY(ILIROK) = YARR(LOCNXT)
            ELSE
               ILDIAG = ILOC(K - 1,K - 1)
               ILOTHR = ILOC(IROW,K - 1)
               POLY(ILIROK) = (POLY(ILDIAG)*DX(IROW)-POLY(ILOTHR)*DX(K-1))/(X(&
                  IROW)-X(K-1))
            ENDIF
         END DO
!
!   Pick off the diagonal element
!
         ILDIAG = ILOC(IROW,IROW)
         EST(IROW) = POLY(ILDIAG)
!
      END DO
!
!   Now the estimate vector is filled in, so obtain the
!   best estimate
!
      DEBEB = ABS((EST(2)-EST(1))/EST(2))

      IBEST = 2
      DO IROW = 3, MXORD
         DEBE = ABS((EST(IROW)-EST(IROW-1))/EST(IROW))

!ps 13/12/2017
!ps If NARR is small, the EST array sometimes contains the NaN values.
!ps In a consequence, DEBE becomes NaN, and in that case the condition
!ps in the CYCLE statement (NaN >= number) is then evaluated to "false".
!ps So, the CYCLE is not performed, what is incorrect.
!ps To avoid this, the CYCLE was replaced with the IF / THEN / END IF
!ps (the condition was reversed).
        !ps IF (DEBE >= DEBEB) CYCLE
        !ps DEBEB = DEBE
        !ps IBEST = IROW
        IF (DEBE < DEBEB) THEN
            DEBEB = DEBE
            IBEST = IROW
        END IF

      END DO
      YVAL = EST(IBEST)
!
      IF (DEBEB > ACCY) WRITE (*, 301) DEBEB, ACCY
!
      RETURN
!
  300 FORMAT('INTERP: Extrapolating, not interpolating.')
  301 FORMAT('INTERP: Accuracy of interpolation (',1P,1D10.3,') is',&
         ' below input criterion (',1D10.3,').')
      RETURN
      CONTAINS


      INTEGER FUNCTION ILOC (IND1, IND2)
      INTEGER, INTENT(IN) :: IND1
      INTEGER, INTENT(IN) :: IND2
      ILOC = (IND1*(IND1 - 1))/2 + IND2
      RETURN
      END FUNCTION ILOC
!
      END SUBROUTINE INTERP
