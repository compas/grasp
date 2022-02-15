!***********************************************************************
!                                                                      *
      FUNCTION SYSTEMMem()
!                                                                      *
!***********************************************************************
!...Written by Chongyang Chen, Fudan University, Shanghai, China       *
!                                                Oct. 2020             *
!***********************************************************************
!...Translated by Gediminas Gaigalas  May 2021
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!MemTotal:       264566412 kB
!MemTotal:       1056570236 kB
      REAL(DOUBLE) :: SYSTEMMem
      INTEGER      :: SUMMem, LENG
      CHARACTER(LEN=128) :: LINE,CMem
!-----------------------------------------------
!
      OPEN(5321,FILE='/proc/meminfo', STATUS='OLD', FORM='FORMATTED')
      READ(5321,'(A)')LINE
      CLOSE(5321)
      LENG = LEN_TRIM(LINE)
      CMem = LINE(10:LENG-3)
      READ(CMem,*)SUMMem  ! KB
!GG      SYSTEMMem = SUMMem / 1.0D6  ! GB
      SYSTEMMem = SUMMem / 1048576  ! GB

      RETURN
      END FUNCTION SYSTEMMem
