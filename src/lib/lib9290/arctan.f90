!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION ARCTAN (ARG1, ARG2)
!-----------------------------------------------
!                                                                      *
!                   -1                                                 *
!       ARCTAN = tan   (ARG1/ARG2),      0 .LE. ARCTAN .LT. 2*\pi      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!   Editted by C. Froese Fischer after translation                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:34   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEF_C
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN) :: ARG1
      REAL(DOUBLE), INTENT(IN) :: ARG2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL, SAVE :: FIRST = .TRUE., INTRIN = .TRUE.
!-----------------------------------------------
!
!
!   Determine whether the FORTRAN intrinsic function ATAN2 always
!   returns a positive value
!
      IF (FIRST) THEN
         ARCTAN = ATAN2(-1.0D00,-1.0D00)
         IF (ARCTAN > 0.0D00) THEN
            INTRIN = .TRUE.
         ELSE
            INTRIN = .FALSE.
         ENDIF
         FIRST = .FALSE.
      ENDIF
!
!   Use the intrinsic function if it passes the above test; otherwise
!   add 2*PI to the negative values returned by the intrinsic function
!
      IF (INTRIN) THEN
         ARCTAN = ATAN2(ARG1,ARG2)
      ELSE
         IF (ARG1 >= 0.0D00) THEN
            ARCTAN = ATAN2(ARG1,ARG2)
         ELSE
            ARCTAN = PI + PI + ATAN2(ARG1,ARG2)
         ENDIF
      ENDIF
!
      RETURN
      END FUNCTION ARCTAN
