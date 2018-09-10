      MODULE perko1_I
      INTERFACE
!
      SUBROUTINE PERKO1(JA,BK,IK,BD,ID)
      USE vast_kind_param, ONLY:  DOUBLE
      INTEGER,      INTENT(IN)                :: JA
      INTEGER,      INTENT(OUT), DIMENSION(7) :: IK, ID
      REAL(DOUBLE), INTENT(OUT), DIMENSION(3) :: BK, BD
      END SUBROUTINE
      END INTERFACE
      END MODULE
