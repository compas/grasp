!***********************************************************************
!
      SUBROUTINE INIESTmpi (NMAX, NCF, NIV, BASIS, HMX, JCOL, IROW) 
!
!       Routine for providing initial estimates from the diagonal 
!     of the matrix. This way was used by Dvdson in atomic structure 
!     calculations. It should be used to obtain estimates when nothing 
!     else is available.
!
!       nmax is typically 1000
 
!  Structure of the input sparse matrix hmx:
!    . It's a 1-d array
!    . Length: 1 to jcol(ncf)
!    . Number of non-zero elements for column j is:
!          jcol(j) - jcol(j-1) + 1
!    . Row index for element hmx(i) is irow(i)
!  Xinghong He  98-10-28
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:36:59   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE memory_man
      use mpi_C

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NMAX 
      INTEGER, INTENT(IN) :: NCF 
      INTEGER  :: NIV 
      INTEGER, INTENT(IN) :: JCOL(0:*) 
      INTEGER, INTENT(IN) :: IROW(*) 
      real(double)  :: BASIS(*) 
      real(double), INTENT(IN) :: HMX(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS, JOFFSPAR, J, JOFFNORM, IR, NFOUND, INFO

      integer, dimension(:), pointer :: iwork,ifail 
      real(double), dimension(:), pointer :: ap, eigval,vec, work 
      real(double) :: ABSTOL

      REAL(KIND(0.0D0)) :: dlamch
!     real dLAMCH
!     external dLAMCH

!-----------------------------------------------

      NS = MIN(NMAX,NCF) 
      CALL ALLOC (AP, (NS*(NS + 1))/2, 'AP', 'INIESTmpi')
      CALL DINIT ((NS*(NS + 1))/2, 0.D0, AP, 1) 
 
!  Expand the sparse form to normal form for upper-right sub-matrix
 
      JOFFSPAR = 0                               ! offset for sparse form 
!     DO J = 1, NS 
      DO j = myid + 1, ns, nprocs
         JOFFNORM = (J*(J - 1))/2                ! offset for normal form 
         DO IR = JOFFSPAR + 1, JCOL(J) 
            AP(IROW(IR)+JOFFNORM) = HMX(IR) 
         END DO 
         JOFFSPAR = JCOL(J) 
      END DO 
 
!  Merge ap from all nodes and then send to all nodes
 
      CALL gdsummpi (ap, (NS*(NS+1))/2)

      CALL ALLOC(eigval, ns,     'EIGVAL', 'INIESTmpi')
      CALL ALLOC (vec,   ns*niv, 'VEC',    'INIESTmpi')
      CALL ALLOC (work,  8*ns,   'WORK',   'INIESTmpi')
!GG      CALL ALLOC (iwork, 8*ns,   'IWORK',  'INIESTmpi')
      CALL ALLOC (iwork, 5*ns,   'IWORK',  'INIESTmpi')
      CALL ALLOC (ifail, ns,     'IFAIL',  'INIESTmpi')
 
      ABSTOL = 2*dLAMCH('S')

      CALL DSPEVX ('Vectors also', 'In a range', 'Upper triangular', NS, AP, &
      -1., -1., 1, NIV, ABSTOL, NFOUND, EIGVAL, VEC, NS, WORK, IWORK, IFAIL, &
         INFO) 
      IERR = -ABS(INFO) 
 
!  Build the Basis.
 
      CALL DINIT (NCF*NIV, 0.D0, BASIS, 1) 
 
!  scatter the vectors
 
      DO J = 1, NIV 
         CALL DCOPY (NS, VEC(NS*(J-1)+1), 1, BASIS(NCF*(J-1)+1), 1) 
      END DO 
 
      CALL DCOPY (NIV, EIGVAL, 1, BASIS(NIV*NCF+1), 1) 
 
       CALL DALLOC (ap,     'AP',     'INIESTmpi')
       CALL DALLOC (eigval, 'EIGVAL', 'INIESTmpi')
       CALL DALLOC (vec,    'VEC',    'INIESTmpi')
       CALL DALLOC (work,   'WORK',   'INIESTmpi')
       CALL DALLOC (iwork,  'IWORK',  'INIESTmpi')
       CALL DALLOC (ifail,  'IFAIL',  'INIESTmpi')
 
      RETURN  
      END SUBROUTINE INIESTmpi
