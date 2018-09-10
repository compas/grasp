      MODULE mpi_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!cjb
!   "INCLUDE mpif.h" is going to be deprecated in future MPI releases
!   and "USE mpi" is recommended
!   if you already have an MPI library with precompiled mpi.mod
!   we recommend to uncomment version "USE mpi" below 
!   and comment out version 'mpif.h'
!   depending on your environment
!   it may be necessary to add $(MPI_INC) = -I/path  in Makefile
!   grasp/src/lib/mpi90/Makefile
!   with 'path' pointing to module 'mpi.mod'
!cjb
!   otherwise uncomment 'mpif.h' and comment out "USE mpi" 
!   do not uncomment both
!cjb

!cjb begin version USE mpi
! uncomment the next 2 lines
!     USE mpi
!     implicit none
!cjb end   version USE mpi

!cjb begin version INCLUDE 'mpif.h'
! uncomment the next 2 lines
      implicit none
      INCLUDE 'mpif.h'
!cjb end   version INCLUDE 'mpif.h'

      character*(MPI_MAX_PROCESSOR_NAME) host
      INTEGER myid, nprocs, lenhost
      INTEGER ierr, ierrtotal
      INTEGER istat(MPI_STATUS_SIZE)

      END MODULE mpi_C
