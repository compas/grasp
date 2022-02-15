!***********************************************************************
!                                                                      *
      subroutine hmoutmpi(myid, nprocs, ncf)
!                                                                      *
!   Routine for printing the Hamiltonian matrix.                       *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
!   Block Version by Xinghong He          Last revision: 30 Jan 1999   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
      USE hmat_C
!GG      USE mpi_C
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: myid
      integer  :: nprocs
      integer, intent(in) :: ncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ibeg, ico, idiag, list, iro
!-----------------------------------------------
!
      ibeg = 1
      do ico = myid + 1, ncf, nprocs
         idiag = iendc(ico)
         do list = ibeg, idiag
            iro = irow(list)
            write (99,*) 'H(',iro,ico,')= ', emt(list)
         enddo
         ibeg = idiag + 1
      enddo
!
      return
      end subroutine hmoutmpi
