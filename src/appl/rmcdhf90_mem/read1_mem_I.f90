      MODULE read1_mem_I
      INTERFACE
!   Written by  G. Gaigalas               Vilnius, September 2021
      SUBROUTINE read1_mem(NFILE, NCONTR_tot2, LAB, NCONTR)
      USE vast_kind_param, ONLY: LONG
      INTEGER(LONG), INTENT(IN) :: NCONTR_tot2
      INTEGER, INTENT(IN)       :: NFILE
      INTEGER, INTENT(OUT)      :: LAB
      INTEGER, INTENT(OUT)      :: NCONTR
      END SUBROUTINE
      END INTERFACE
      END MODULE
