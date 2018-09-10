      MODULE cle0sm_I
      INTERFACE
!
      SUBROUTINE CLE0SM(Q,QM,S,C,CM,A)
      USE vast_kind_param, ONLY:  DOUBLE
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, S, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
      END SUBROUTINE
      END INTERFACE
      END MODULE
