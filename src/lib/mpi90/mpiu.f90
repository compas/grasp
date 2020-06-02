! Collection of FORTRAN subroutines for general purpose MPI operations.
! Mostly IMPLICIT NONE
! Only INCLUDE 'mpif.h'
!***********************************************************************
      subroutine startmpi (myid, nprocs, host, lenhost)
! This subroutine starts MPI to get the process id (myid), the number
! of processes (nprocs) and the status of enrollment (ierr)

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: myid,  nprocs, lenhost
      CHARACTER(LEN=*), INTENT(IN) :: host
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      include 'mpif.h'
      integer ierr

      call MPI_Init (ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, nprocs, ierr)
      call MPI_Get_processor_name (host, lenhost, ierr)
      if (ierr .ne. MPI_SUCCESS) print *, ' ierr=', ierr, ' id=', myid

      return
      end

!***********************************************************************
      subroutine startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, progname)
! Calls startmpi to get mpi environment;
! Calls cpath to get various paths;
! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: myid,  nprocs, lenhost
      INTEGER, INTENT(OUT)  :: ncount1
      CHARACTER(LEN=*), INTENT(IN) :: host, progname
      CHARACTER(LEN=*), INTENT(OUT) :: startdir, permdir, tmpdir
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      character idstring*3
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond
      integer lenstart, lenperm, lentmp, ncount_rate, ncount_max
! For printing
      CHARACTER(LEN=80) :: msg

!=======================================================================
!  Get processor info: myid, nprocs, host name; and print
!=======================================================================

      CALL startmpi (myid, nprocs, host, lenhost)
      WRITE (idstring, '(I3.3)') myid
      IF (myid .EQ. 0) THEN
         print *, '===================================================='
         print *, '       ', progname, ': Execution Begins ...'
         print *, '===================================================='
         print *,        'Participating nodes:'
       ENDIF
      msg = '  Host: ' // host(1:lenhost) // '    ID: ' // idstring
      CALL mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Get date, time, zone and print
!=======================================================================

      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)
      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Date and Time:'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' // &
            '  Date: ' // chdate // &
            '  Time: ' // chtime // &
            '  Zone: ' // chzone
      CALL mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Set up local working dir and go there
!     tmpdir  - local working dir of the node. mcpXXX files are there
!     permdir - for I/O specific to node-0.
!=======================================================================

      CALL cpath (startdir, permdir, tmpdir)

      lenstart = LEN_TRIM (startdir)
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Start Dir:'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' // startdir(1:lenstart)
      CALL mpix_printmsg (msg, myid, nprocs)

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Serial I/O Dir (node-0 only):'
         PRINT *, '  ' // host(1:lenhost) // ': ' // permdir(1:lenperm)
      ENDIF

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Work Dir (Parallel I/O):'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' // tmpdir(1:lentmp)
      CALL mpix_printmsg (msg, myid, nprocs)

!=======================================================================
!  Start timing - Record the wall clock
!=======================================================================

      CALL SYSTEM_CLOCK (ncount1, ncount_rate, ncount_max)
      return
      end

!***********************************************************************
      subroutine stopmpi (what, myid)

! This subroutine stops MPI, but before doing so it issues an error
! message to node-0 about what went wrong and where it happened.
!   what - string, type of error (most likely a subroutine name)
!   myid - The id of the PE where error happened

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER, INTENT(IN)  ::  myid
      CHARACTER(LEN=*), INTENT(IN) :: WHAT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ierr

      print *, 'mpi stopped by node-', myid, ' from ', what
      call MPI_Barrier (MPI_COMM_WORLD,ierr)
      call MPI_Finalize (ierr)

      stop
      end

!***********************************************************************
      subroutine stopmpi2 (myid, nprocs, host, lenhost, ncount1, &
                            progname)

! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER, INTENT(IN)  ::  myid,nprocs,lenhost
      integer              ::  ncount1
      character*(*), INTENT(IN) ::  host, progname
! Things for timing
      INTEGER   ncount2, ncount_rate, ncount_max, nseconds
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

! For printing
      character str2nds*8, msg*80

      integer ierr

!=======================================================================
!  Get processor info: myid, nprocs, host name; and print
!=======================================================================

      if (myid .eq. 0) then
         print *, '===================================================='
         print *, '       ', progname, ': Execution Finished ...'
         print *, '===================================================='
         print *,        'Wall time:'
       endif

      call system_clock (ncount2, ncount_rate, ncount_max)
      ncount2 = ncount2 - ncount1
      nseconds = ncount2 / ncount_rate
      write (str2nds, '(i8)') nseconds
      msg = str2nds // ' seconds on ' // host(1:lenhost)
      call mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) then
         print *
         print *, 'Finish Date and Time:'
      endif

      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  ' // host(1:lenhost) // ': ' // &
            '  Date: ' // chdate // &
            '  Time: ' // chtime // &
            '  Zone: ' // chzone
      CALL mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) print *

      call MPI_Barrier (MPI_COMM_WORLD,ierr)
      call stopmpi (progname // ': Execution complete.', myid)
      return
      end

!***********************************************************************
      subroutine mpix_printmsg (msg, myid, nprocs)

! Displays on node-0's screen info from all nodes including node 0.

!***********************************************************************
!  Modified by Charlotte F. Fischer   10/10/2017

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER, INTENT(IN)  ::  myid, nprocs
      CHARACTER(LEN=*), INTENT(OUT) :: msg

!-----------------------------------------------
!   L o c a l  V a r i a b l e s
!-----------------------------------------------
      INTEGER :: inID, istat(MPI_STATUS_SIZE), ierr, msgLength

      msgLength = len_trim (msg)

      if (myid .ne. 0) then
         call MPI_Send (msgLength, 1, MPI_INTEGER, 0, myid, &
                        MPI_COMM_WORLD, ierr)   ! Send nsgLength
         call MPI_Send (msg, msgLength, MPI_CHARACTER, 0, myid+nprocs, &
                        MPI_COMM_WORLD, ierr)   ! Send msg
      else
         print *, msg(1:msgLength)      ! msg from node 0 itself
         do inID = 1, nprocs - 1
            call MPI_Recv (msgLength, 1, MPI_INTEGER, inID, &
                           inID, MPI_COMM_WORLD, istat, ierr)
            call MPI_Recv (msg, msgLength, MPI_CHARACTER, inID, &
                           inID+nprocs, MPI_COMM_WORLD, istat, ierr)
            print *, msg(1:msgLength)
         enddo
      endif

      return
      end

!***********************************************************************
      subroutine mpix_chkpt (myid, what)
!     To set a check-point in mpi program
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER, INTENT(IN) :: myid
      CHARACTER(LEN=*), INTENT(IN) :: what
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ierr
      CHARACTER(LEN=1) :: CHTMP
!-----------------------------------------------
!
      if (myid .eq. 0) then
         print *, what, 'hit any key to continue'
         read (*, '(a)') chtmp
      endif

      call mpi_barrier (MPI_COMM_WORLD, ierr)

      return
      end

!***********************************************************************
      subroutine mpix_bytes (n, newType, ierr)

!  Constructs new mpi data type newType of n-bytes long
!***********************************************************************
!  Modified by Charlotte F. Fischer  10/10/2017
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      IMPLICIT NONE
      include 'mpif.h'
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(out) :: newtype, ierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ier0

      call MPI_Type_contiguous (n, MPI_BYTE, newType, ier0)
      call MPI_Type_commit (newType, ierr)

      ierr = ierr + ier0

      return
      end

!***********************************************************************
      SUBROUTINE gdsummpi (x, n)

!     Sum x over all nodes and in y, then copy y to  x
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN)                         :: n
      REAL(DOUBLE), DIMENSION(1:n), INTENT(INOUT) ::   x

      INTEGER                                     :: ierr
      REAL(DOUBLE), DIMENSION(1:n)                :: y

      CALL dinit (n, 0.d0, y, 1)
      CALL MPI_Allreduce (x, y, n, MPI_DOUBLE_PRECISION, &
                                   MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (n, y, 1, x, 1)   ! copy y to x

      RETURN
      END

!***********************************************************************
      SUBROUTINE gisummpi (ix, n)
!     Sum x over all nodes and in iy, then copy iy to  ix
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN)                         :: n
      INTEGER, DIMENSION(1:n), INTENT(INOUT) ::  ix

      INTEGER                                     :: ierr, i
      INTEGER, DIMENSION(1:n)                ::  iy

      iy = 0

      CALL MPI_Allreduce (ix, iy, n, MPI_INTEGER, MPI_SUM, &
                                     MPI_COMM_WORLD, ierr)
      CALL icopy (n, iy, 1, ix, 1)

      RETURN
      END

!***********************************************************************
!cychen, 2020/05
      SUBROUTINE ddotmpi (n, x, y, ddot)
!     mpi inner-dot, modification from blas: ddot
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mpi_C
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER n, i, nbar
      REAL*8  x(n), y(n), ddot, dsum

      ddot=0.0d0
      nbar = n/nprocs

      do i=myid*nbar+1, (myid+1)*nbar
         ddot=ddot+x(i)*y(i)
      enddo

      CALL MPI_Allreduce (ddot, dsum, 1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
      do i=nbar*nprocs+1, n
         dsum=dsum+x(i)*y(i)
      enddo
      ddot=dsum

      RETURN
      END

!***********************************************************************
!cychen, 2020/05
      subroutine dgemvmpi ( trans, m, n, alpha, a, lda, x, incx, &
                        BETA, Y, INCY )
!     mpi version for blas: dgemv
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mpi_C
!-----------------------------------------------
      IMPLICIT NONE
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )

      INTEGER  nbar

      real*8   tm(m), tn

!cychen:
!DGEMV performs one of the matrix-vector operations by using mpi:
! y :=A*x, or y:=A'*x
!then, alpha and beta must be input as 1.0d0 and 0.0d0, respectively
!Also, for simplicity, incx and incy are modified to be both 1.

!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
              .NOT.LSAME( TRANS, 'T' ).AND.  &
              .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
!cychen      ELSE IF( INCX.EQ.0 )THEN
      ELSE IF( INCX.NE.1 )THEN
         INFO = 8
         if (myid.eq.0) print *, "INCX=",INCX
!cychen      ELSE IF( INCY.EQ.0 )THEN
      ELSE IF( INCY.NE.1 )THEN
         INFO = 11
         if (myid.eq.0) print *, "INCY=",INCY
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMVMPI ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
         ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
        RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
!cychen      IF( BETA.NE.ONE )THEN
!cychen         IF( INCY.EQ.1 )THEN
!cychen            IF( BETA.EQ.ZERO )THEN
!cychen               DO 10, I = 1, LENY
!cychen                  Y( I ) = ZERO
!cychen   10          CONTINUE
!cychen            ELSE
!cychen               DO 20, I = 1, LENY
!cychen                  Y( I ) = BETA*Y( I )
!cychen   20          CONTINUE
!cychen            END IF
!cychen         ELSE
!cychen            IY = KY
!cychen            IF( BETA.EQ.ZERO )THEN
!cychen               DO 30, I = 1, LENY
!cychen                  Y( IY ) = ZERO
!cychen                  IY      = IY   + INCY
!cychen   30          CONTINUE
!cychen            ELSE
!cychen               DO 40, I = 1, LENY
!cychen                  Y( IY ) = BETA*Y( IY )
!cychen                  IY      = IY           + INCY
!cychen   40          CONTINUE
!cychen            END IF
!cychen         END IF
!cychen      END IF
!cychen      IF( ALPHA.EQ.ZERO )
!cychen     $   RETURN

      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
!         JX = KX
         IF( INCY.EQ.1 )THEN
            CALL DINIT(M, 0.d0, tm, 1)
            nbar = M / nprocs
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                do i = myid * nbar + 1, (myid + 1) * nbar
                   tm(i)=tm(i) + x(j) * A(i,j)
                enddo
               END IF
   60       CONTINUE
            CALL MPI_Allreduce (tm, y, m, MPI_DOUBLE_PRECISION, &
                                  MPI_SUM, MPI_COMM_WORLD, ierr)
            do j = 1, N
               do i = nbar*nprocs + 1, M
                  y(i) = y(i) + x(j) * A(i,j)
               enddo
            enddo

         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            nbar = M / nprocs
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = myid * nbar + 1, (myid + 1) * nbar
                  TEMP = TEMP + A( I, J ) * X( I )
   90          CONTINUE
               CALL MPI_Allreduce (TEMP, tn, 1, MPI_DOUBLE_PRECISION,&
                                  MPI_SUM, MPI_COMM_WORLD, ierr)
               do i = nbar * nprocs + 1, M
                  tn = tn + A( I, J ) * X( I )
               enddo
               Y( J ) = tn
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END

!***********************************************************************
!cychen, 2020/05
!performs y := da*x + y
!For simplicity, modification only for incx=incy=1.
      subroutine daxpympi(n,da,dx,incx,dy,incy)
!     mpi version of blas: daxpy (Tests show that it does not save time)
!***********************************************************************
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE mpi_C
!-----------------------------------------------

      IMPLICIT NONE

      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n

      INTEGER  nbar
      real*8   tm(n)

!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

!     mpi code for both increments equal to 1
   20 CALL DINIT(n, 0.d0, tm, 1)
      nbar = n / nprocs
      do i = myid*nbar+1, (myid+1)*nbar
         tm(i) = dy(i) + da * dx(i)
      enddo
      CALL MPI_Allreduce (tm, dy, nprocs*nbar, MPI_DOUBLE_PRECISION,&
                        MPI_SUM, MPI_COMM_WORLD, ierr)
      do i = nprocs*nbar+1, n
         dy(i) = dy(i) + da*dx(i)
      enddo

!        m = mod(n,4)
!      if( m .eq. 0 ) go to 40
!      do 30 i = 1,m
!        dy(i) = dy(i) + da*dx(i)
!   30 continue
!      if( n .lt. 4 ) return
!   40 mp1 = m + 1
!      do 50 i = mp1,n,4
!        dy(i) = dy(i) + da*dx(i)
!        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
!        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
!        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
!   50 continue
      return
      end
