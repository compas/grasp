      MODULE read2_mem_I
      INTERFACE
!   Written by  G. Gaigalas               Vilnius, September 2021
      SUBROUTINE read2_mem(NFILE,NCONTR_tot,NCONTR,ICLMNDUM)
      USE vast_kind_param, ONLY: LONG
      INTEGER(LONG), INTENT(IN) :: NCONTR_tot
      INTEGER, INTENT(IN)       :: NFILE
      INTEGER, INTENT(OUT)      :: NCONTR
      INTEGER, INTENT(OUT)      :: ICLMNDUM
      END SUBROUTINE
      END INTERFACE
      END MODULE
