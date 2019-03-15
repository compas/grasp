!***********************************************************************
!                                                                      *
      SUBROUTINE allocCheck(n, IREZ)
!                                                                      *
!   Written by G. Gaigalas            Last revision: 27 October 2017   *
!                                                                      *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: IREZ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), POINTER, DIMENSION(:) :: pr
      INTEGER, POINTER, DIMENSION(:) :: pi
      INTEGER :: error
!-----------------------------------------------
      IREZ = 1
      allocate(pr(n),STAT=error)
      IF (error /= 0) then
         IREZ = 0
         RETURN
      END IF
      allocate(pi(3*n),STAT=error)
      IF (error /= 0) then
         IREZ = 0
         CALL DALLOC (PR,'PR','allocCheck')
         RETURN
      END IF
      CALL DALLOC (PR,'PR','allocCheck')
      CALL DALLOC (PI,'PI','allocCheck')
      RETURN
      END SUBROUTINE allocCheck
