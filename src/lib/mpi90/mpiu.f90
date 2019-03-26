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
