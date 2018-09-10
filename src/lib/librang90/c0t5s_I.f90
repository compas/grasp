      MODULE c0t5s_I
      INTERFACE
!
      SUBROUTINE C0T5S(Q, QM, SM, C, CM, A)
      USE vast_kind_param, ONLY:  DOUBLE
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, SM, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
      END SUBROUTINE
      END INTERFACE
      END MODULE
