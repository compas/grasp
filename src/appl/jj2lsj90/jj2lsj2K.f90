!
!***********************************************************************
!                                                                      *
      PROGRAM jj2lsj2K
!                                                                      *
!     This MAIN program controls the transformation of atomic states,  *
!     which are given in a jj-coupled CSF basis, into an LS-coupled    *
!     basis. The program requires a jj-coupled basis in standard order *
!     where, if both subshells | n j = l-1/2> and | n j = l+1/2>       *
!     of a given shell (nl) occurs, they always follow successively    *
!     in this order. The LS-coupled basis, moreover, is given in       *
!     the same sequence of shells.                                     *
!                                                                      *
!     All LS-jj transformation coefficients are precalculated and      *
!     'stored' in the modules rabs_lsj_data_1, rabs_lsj_data_2 and     *
!     rabs_lsj_data_3.                                                 *
!                                                                      *
!     Calls:  FACTT, SETISO, JJ2LSJ, starttime, stoptime.              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  May 2011   *
!     VILNIUS                                               May 2017   *
!                                                                      *
!     Modified by G. Gaigalas and C. Cychen                     2021   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE jj2lsj_code
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer           :: ncount1
!-----------------------------------------------
      call starttime (ncount1, 'jj2lsj')
      print *, " "
      print *, "jj2lsj: Transformation of ASFs from a jj-coupled CSF basis"
      print *, "        into an LS-coupled CSF basis  (Fortran 95 version)"
      print *, "        (C) Copyright by   G. Gaigalas and Ch. F. Fischer,"
      print *, "        (2021)."
      print *, "        Input files: name.c, name.(c)m"
      print *, "        Ouput files: name.lsj.lbl,"
      print *, "         (optional)  name.lsj.c, name.lsj.j,"
      print *, "                     name.uni.lsj.lbl, name.uni.lsj.sum"
      print *, " "
!
!  Set up the table of logarithms of factorials
      call factt
!
      CALL setiso('isodata')
      CALL jj2lsj
!
      call stoptime (ncount1, 'jj2lsj')
      stop
      end program jj2lsj2K
