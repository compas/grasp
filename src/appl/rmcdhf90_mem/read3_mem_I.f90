      MODULE read3_mem_I
      INTERFACE
!   Written by  G. Gaigalas               Vilnius, September 2021
      SUBROUTINE read3_mem(NFILE,NCONTR_tot,NCONTR)
      USE vast_kind_param, ONLY: LONG
      INTEGER(LONG), INTENT(IN) :: NCONTR_tot
      INTEGER, INTENT(IN)       :: NFILE
      INTEGER, INTENT(OUT)      :: NCONTR
      END SUBROUTINE
      END INTERFACE
      END MODULE
