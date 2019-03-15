!************************************************************************
!
      SUBROUTINE INIESTdm (nmax, ncf, NIV, BASIS, hmx)

!     Routine for providing initial estimates from the diagonal
!     of the matrix. This way was used by Dvdson in atomic structure
!     calculations. It should be used to obtain estimates when nothing
!     else is available.
!
!   Block version by Xinghong He            Last revision: 18 Jun 1998
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE eigv_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dspevx_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN):: nmax, ncf, niv
      REAL(DOUBLE), DIMENSION(*) :: basis, hmx
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      INTEGER :: info, ns, i, nfound, j, ihmx, iap
      REAL(DOUBLE), DIMENSION(:), pointer :: ap, eigval, vec, work
      INTEGER, DIMENSION(:), pointer :: iwork, ifail
      REAL(DOUBLE) :: ABSTOL
      REAL :: SLAMCH
!-----------------------------------------------------------------------
      NS = MIN (nmax, ncf)

      CALL alloc (ap, NS*(NS+1)/2, 'AP', 'INIESTDM')
      CALL dinit (NS*(NS+1)/2, 0.d0, ap, 1)

!  Get the upper left sub-matrix
!  ap is global, hmx is local !!!

      ihmx = 0
      DO i = myid + 1, ns, nprocs
         iap = i*(i-1)/2
         DO j = 1, i
            ap(j+iap) = hmx(j+ihmx)
         ENDDO
         ihmx = ihmx + i
      ENDDO

      CALL alloc (eigval,   NS, 'EIGVAL', 'INIESTDM')
      CALL alloc (vec,   NS*NIV, 'VEC', 'INIESTDM')
      CALL alloc (work,  8*NS,  'WORK', 'INIESTDM')
      CALL alloc (iwork, 5*NS,  'IWORK', 'INIESTDM')
      CALL alloc (ifail,     NS,  'IFAIL', 'INIESTDM')

!      CALL DSPEVX ('Vectors also','In a range','Upper triangular',     &
!     &          NS,AP,-1.0D0,-1.0D0,1,NIV,0.d0,                        &
!     &          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      ABSTOL = 2*SLAMCH('S')
      CALL DSPEVX ('V','I','U',                                        &
                NS,AP,-1.d0,-1.d0,1,NIV,ABSTOL,                            &
                NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      IERR = -ABS (INFO)

!  Build the Basis.

      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)

!  scatter the vectors

      DO J = 1, NIV
            CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO

      CALL dcopy (NIV, EIGVAL, 1, BASIS(NIV*ncf+1), 1)

      CALL dalloc (ap, 'AP', 'INIESTDM')
      CALL dalloc (eigval, 'EIGVAL', 'INIESTDM')
      CALL dalloc (vec, 'VEC', 'INIESTDM')
      CALL dalloc (work, 'WORK', 'INIESTDM')
      CALL dalloc (iwork, 'IWORK', 'INIESTDM')
      CALL dalloc (ifail, 'IFAIL', 'INIESTDM')

      RETURN
      END
