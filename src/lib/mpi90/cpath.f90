!***********************************************************************
      subroutine cpath (startdir, permdir, tmpdir)
!
!      startdir - path where the current node started from.
!      permdir  - path where node-0 performs serial i/o.
!      tmpdir   - path where the current node goes to.
!
!   This version reads (by node-0) the paths from a disk file under the
!   starting directory of the node-0, determine the length and do
!   sending/receiving. Only if the paths defined here do not exist will
!   C functions be called to create them.
!
!   Xinghong He  98-10-30
!
!   updated for a case when tmpdir did not exist prior to execution, and
!   has/have to be created with call sys_mkdir
!                                               Jacek Bieron 2017-10-31
!
!   cpath makes a sequence of attempts to create tmpdir directories
!   at local disks of all nodes, in the following order:
!   (1) 'disks' file
!   (2) env variable MPI_TMP
!   (3) directory /scratch/$USER
!     if file 'disks' exists, cpath reads it
!     if file 'disks' does not exist, cpath uses env variable MPI_TMP
!     if env variable MPI_TMP is not set, cpath defaults to /scratch/$USER
!   if error occurs, cpath tries twice more
!   if error persists, cpath attempts to cleanly terminate MPI_COMM_WORLD
!   by a call to stopmpi with error message or by invoking MPI_ABORT
!
!cjb                                            Jacek Bieron 2018 May 10
!   numbering of comments:
!   even numbers = 22 steps
!   odd numbers = errors and warnings
!
!   function 'system' and execute_command_line commented out
!cjb                                            Jacek Bieron 2018 June 18
!
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!     COMMON /mpi/ myid, nprocs, ierr
      use mpi_C
!     include 'mpif.h'
!     integer  myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
      character(len=*), intent(out):: startdir, permdir, tmpdir
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  ierr_all_nodes, ierr2m
!     ierr is used by cpath calls to sys_*dir routines
!     ierr2m is used by MPI routines and by 'system' routines
      !...locals - be careful when changing lenidstring !
      integer, parameter:: lendisk0 = 128, lenidstring = 3
      character(len=lendisk0) disk
      character(len=lendisk0) mpi_tmp
      character(len=lenidstring) idstring
      integer   lendisk, i, len_cwd
      character*(64) en,ev,uname
      integer lstring
      integer asleep, itry, iteery
      logical ydisks
!     integer system
!    integer lengthsleep, statusleep

!=======================================================================
! step 0 get working directory
!=======================================================================
      call sys_getwd (startdir,len_cwd)
!     print *, ' in cpath myid = ', myid, ' startdir = ', startdir

!=======================================================================
!  Open disks file, read paths and send/receive them. Each node will have
!  its preliminary path stored in variable disk. In addition, node-0 will
!  have the current working directory stored in permdir
!=======================================================================

!=======================================================================
      if (myid .eq. 0) then
         inquire(FILE='disks',exist=ydisks)
         if(ydisks) then
!=======================================================================
! step 02 read disks file and MPI_Send
!=======================================================================
            open (unit=1001, file='disks', status='old')
            !...paths for serial i/o, node-0 only
            read (1001,*) permdir  ! paths for serial i/o, node-0 only
            read (1001,*) tmpdir   ! temporary for local disk of node-0
            !...paths for slaves, read and send;
            do i = 1, nprocs - 1
               read (1001,*) disk
               call MPI_Send (disk, lendisk0, MPI_CHARACTER, i, i, &
                              MPI_COMM_WORLD, ierr2m)
            enddo
            disk = tmpdir           ! local disk of node-0
            close (1001)
         else
            permdir = startdir
!
! check whether user had set env variable MPI_TMP
! examples:
! csh: setenv MPI_TMP "${HOME}/grasprun/"
! bash: export MPI_TMP = "/scratch/${USER}"
!
!=======================================================================
! step 04 get env variable MPI_TMP
!=======================================================================
             en = "MPI_TMP";
!           call getenv(en,ev)
             call get_environment_variable(en,ev,lstring,ierr)
           if (ierr .ne. 0) then

! error041 failed to get env variable MPI_TMP
      print *, ' warning041 cpath failed at getenv MPI_TMP, myid = ', myid
      print *, ' warning041 cpath will now default to /scratch/USER '
!=======================================================================
! step 06 default to /scratch/USER
! if env variable MPI_TMP is not set, cpath defaults to /scratch/$USER
!=======================================================================
             en = "USER";
             call get_environment_variable(en,ev,lstring,ierr)
             if (ierr .ne. 0) then
! error061 failed to get env variable USER = MPI_ABORT
      print *, ' error061 cpath failed at getenv USER, myid = ', myid
      print *, ' error061 cpath fatal error, calling MPI_ABORT'
               call MPI_ABORT(MPI_COMM_WORLD, ierr, ierr2m)
             else
!=======================================================================
! else back to step 06 default to /scratch/USER
!=======================================================================
               lstring = len_trim(ev);
               mpi_tmp = trim(ev);
               mpi_tmp = '/scratch/'//mpi_tmp(1:lstring)
               lstring = len_trim(mpi_tmp);
             endif
           else
!=======================================================================
! else back to step 04 use env variable MPI_TMP
!=======================================================================
             lstring = len_trim(ev);
             mpi_tmp = trim(ev);
           endif
!
! MPI_TMP/uname
! version MPI_TMP/uname --- commented out
! set tmpdir = MPI_TMP/uname
!            en = "USER";
!            call getenv(en,ev);
!            lstring2 = len_trim(ev);
!            uname = trim(ev);
!           tmpdir = mpi_tmp(1:lstring1)//uname(1:lstring2)

!=======================================================================
! step 08 MPI_Send/MPI_Recv tmpdir
!=======================================================================
            tmpdir = mpi_tmp(1:lstring)
            disk = tmpdir(1:lstring)
            do i = 1, nprocs - 1
               call MPI_Send (disk, lendisk0, MPI_CHARACTER, i, i, &
                              MPI_COMM_WORLD, ierr2m)
            enddo
         end if
      else  ! node0 has sent
         !...slaves receive their local dirs
         if (nprocs .gt. 1) then
           call MPI_Recv (disk, lendisk0, MPI_CHARACTER, 0, myid, &
                          MPI_COMM_WORLD, istat, ierr2m)
         endif
      end if  ! node0 has sent & slaves have received
      lendisk = len_trim (disk)
!
!=======================================================================
!  Each node goes to local disk. Create if did not exist
!=======================================================================

!=======================================================================
!  this section should be (finished and) implemented
!  if cpath does not work properly on a shared memory system due to
!  too many simultaneous disks requests from sys_chdir or sys_mkdir
!
!  sequentially call sys_mkdir (to avoid simultaneous sys_mkdir)
!=======================================================================
!
!=======================================================================
! step 14 mkdir tmpdir
!=======================================================================
!
      call sys_chdir (disk, lendisk, ierr)
!     print *, ' in cpath myid = ', myid, ' disk = ', disk
      if (ierr.ne.0) call sys_mkdir (disk, lendisk, ierr)
      if (ierr.ne.0) then
! error141 failed to mkdir tmpdir
! step14 try mkdir twice more
        do itry = 1, 2
! get some sleep -- different loop times depending on myid
          do iteery = 1, myid+1
            asleep = 2
            asleep = asleep * int(sqrt(real(3*iteery)))
            asleep = asleep / int(sqrt(real(5*iteery)))
          enddo
! get more sleep -- different sleep times (1-9 sec) depending on myid
!         asleep = int(1 + 4 * itry * real(myid) / real(nprocs))
!         ierr2m = system ('sleep ' // achar(asleep+48))
! use EXECUTE_COMMAND_LINE if function 'system' is not supported
!       call execute_command_line ('sleep ' // achar(asleep+48), &
!                                   exitstat=ierr2m)
!       print *, ' error141 cpath achar(asleep+48) = ', achar(asleep+48)
!       print *, ' step14 cpath sleeps ', achar(asleep+48), ' seconds,',&
!                ' myid = ', myid

          call sys_mkdir (disk, lendisk, ierr)
          if (ierr.eq.0) exit
        enddo
!      if (ierr.ne.0) call exit(1)
        if (ierr .ne. 0) then
!
! error141 failed to mkdir tmpdir
        print *, ' error141 cpath failed at sys_mkdir(disk), myid = ',&
                                                             myid
!        stop
!        stop
!        CALL stopmpi (' error in cpath ', myid)
         goto 999
        endif
      endif
!
!=======================================================================
! step 16 cd tmpdir
!=======================================================================
!
      call sys_chdir (disk, lendisk, ierr)
      if (ierr .ne. 0) then
! error161 failed to chdir tmpdir
! step16 try chdir twice more
        do itry = 1, 2
! get some sleep -- different loop times depending on myid
          do iteery = 1, myid+1
            asleep = 2
            asleep = asleep * int(sqrt(real(3*iteery)))
            asleep = asleep / int(sqrt(real(5*iteery)))
          enddo

! get more sleep -- different sleep times (1-9 sec) depending on myid
          asleep = int(1 + 4 * itry * real(myid) / real(nprocs))
!         ierr2m = system ('sleep ' // achar(asleep+48))
! use EXECUTE_COMMAND_LINE if function 'system' is not supported
!       call execute_command_line ('sleep ' // achar(asleep+48), &
!                                   exitstat=ierr2m)
!       print *, ' error161 cpath achar(asleep+48) = ', achar(asleep+48)
!       print *, ' step16 cpath sleeps ', achar(asleep+48), ' seconds,',&
!                ' myid = ', myid
          call sys_chdir (disk, lendisk, ierr)
          if (ierr.eq.0) exit
        enddo
!      if (ierr.ne.0) call exit(1)
        if (ierr .ne. 0) then
!
! error161 failed to chdir tmpdir
       print *, ' error161 cpath failed at sys_chdir(disk); myid = ',&
                                                            myid
!         stop
!         stop
!         CALL stopmpi (' error in cpath ', myid)
          goto 999
        endif
      endif
! chdir ok
!     print *, ' in cpath myid = ', myid, ' local disk = ', disk

!=======================================================================
! step 18 cd tmpdir/myid
!  go to sub-dir defined by myid. Create if did not exist
!=======================================================================
!
      write (idstring,'(I3.3)') myid
      call sys_chdir (idstring, lenidstring, ierr)
      if (ierr .ne. 0) then
         call sys_mkdir (idstring, lenidstring, ierr)
        if (ierr .ne. 0) then
!
! error181 failed to mkdir tmpdir/myid
          print *, ' error181 cpath failed to make sub-dir ' // idstring
!         stop
!         stop
          goto 999
        endif
        call sys_chdir (idstring, lenidstring, ierr)
        if (ierr .ne. 0) then
!
! error183 failed to chdir tmpdir/myid
         print *, ' error183 cpath failed to go to sub-dir ' // idstring
!         stop
!         stop
          goto 999
        endif
      endif
  999 continue

!=======================================================================
! step 999 check if all nodes succeeded
! Handle error cases. If not succeeded, invoke stopmpi
!=======================================================================
!
!  collect ierr -> ierr_all_nodes and check if ierr_all_nodes == 0
!       print *, ' step20 cpath ierr = ', ierr, ' myid = ', myid
!
      ierr_all_nodes = ierr
      call MPI_Barrier(MPI_COMM_WORLD, ierr2m)
      CALL MPI_ALLREDUCE (MPI_IN_PLACE,ierr_all_nodes,1, &
                          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr2m)
! step20 ierr_all_nodes
!       if (myid .eq. 0) then
!        print *, ' step20 cpath ierr_all_nodes = ', ierr_all_nodes, &
!                 ' myid = ', myid
!       endif
!
!    CALL MPI_Reduce (ierr, ierr_all_nodes, 1, MPI_INTEGER, MPI_SUM, 0,
!    &         MPI_COMM_WORLD, ierr)
!     CALL MPI_Bcast (ierr_all_nodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
! error201 non zero ierr
      if (ierr_all_nodes .ne. 0) then
        print *, ' error201 cpath ierr = ', ierr, ' myid = ', myid
!
        if (myid .eq. 0) then
          print *, ' error201 cpath terminated at 999, myid = ', myid
        endif
!       stop
!       stop
        CALL stopmpi (' error in cpath ', myid)
      endif

!=======================================================================
! step 2222 all nodes succeeded
! cpath succeeded
!=======================================================================

      tmpdir = disk
!     print *,' step2222 cpath succeeded myid = ',myid, ' tmpdir = ',&
!                                                       tmpdir
      return
      end
