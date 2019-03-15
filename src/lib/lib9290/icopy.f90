      subroutine icopy(n, ix, incx, iy, incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!*************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:29   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer  :: incx
      integer  :: incy
      INTEGER, DIMENSION(*), INTENT(IN) :: IX
      INTEGER, DIMENSION(*), INTENT(OUT) :: IY
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, mp1
!-----------------------------------------------
!
!
      if (n <= 0) return
!
!        code for unequal increments or equal increments
!          not equal to 1
!
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
      m = mod(n,7)
      if (m /= 0) then
         iy(:m) = ix(:m)
         if (n < 7) return
      endif
      mp1 = m + 1
      iy(mp1:((n-mp1+7)/7)*7-1+mp1) = ix(mp1:((n-mp1+7)/7)*7-1+mp1)
      return
      end subroutine icopy
