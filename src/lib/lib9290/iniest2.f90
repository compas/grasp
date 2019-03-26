!***********************************************************************
!
      SUBROUTINE INIEST2(NMAX, NCF, NIV, BASIS, HMX, JCOL, IROW)
!
!  Serial version of iniestmpi.
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
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NMAX
      INTEGER, INTENT(IN) :: NCF
      INTEGER  :: NIV
      INTEGER, DIMENSION(0:*), INTENT(IN) :: JCOL
      INTEGER, DIMENSION(*), INTENT(IN) :: IROW
      REAL(DOUBLE), DIMENSION(*) :: BASIS
      REAL(DOUBLE), DIMENSION(*), INTENT(IN) :: HMX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS, JOFFSPAR, J, JOFFNORM, IR, NFOUND, INFO, IERR

      integer, dimension(:), pointer :: iwork,ifail
      real(double), dimension(:), pointer :: ap, eigval,vec, work

!-----------------------------------------------

      NS = MIN(NMAX,NCF)
      CALL ALLOC (AP, (NS*(NS + 1))/2, 'AP', 'INIEST2')
      CALL DINIT ((NS*(NS + 1))/2, 0.D0, AP, 1)

!  Expand the sparse form to normal form for upper-right sub-matrix

      JOFFSPAR = 0                               ! offset for sparse form
      DO J = 1, NS
         JOFFNORM = (J*(J - 1))/2                ! offset for normal form
         DO IR = JOFFSPAR + 1, JCOL(J)
            AP(IROW(IR)+JOFFNORM) = HMX(IR)
         END DO
         JOFFSPAR = JCOL(J)
      END DO

!  Merge ap from all nodes and then send to all nodes

      CALL ALLOC(eigval, ns, 'EIGVAL', 'INIEST2')
      CALL ALLOC (vec, ns*niv, 'VEC', 'INIEST2')
      CALL ALLOC (work, 8*ns, 'WORK', 'INIEST2')
      CALL ALLOC (iwork, 8*ns, 'IWORK', 'INIEST2')
      CALL ALLOC (ifail, ns, 'IFAIL', 'INIEST2')

      CALL DSPEVX ('Vectors also', 'In a range', 'Upper triangular', NS, AP, &
         -1., -1., 1, NIV, 0.D0, NFOUND, EIGVAL, VEC, NS, WORK, IWORK, IFAIL, &
         INFO)
      IERR = -ABS(INFO)

!  Build the Basis.

      CALL DINIT (NCF*NIV, 0.D0, BASIS, 1)

!  scatter the vectors

      DO J = 1, NIV
         CALL DCOPY (NS, VEC(NS*(J-1)+1), 1, BASIS(NCF*(J-1)+1), 1)
      END DO

      CALL DCOPY (NIV, EIGVAL, 1, BASIS(NIV*NCF+1), 1)

      !deallocate(ap)
      !deallocate(eigval)
      !deallocate(vec)
      !deallocate(work)
      !deallocate(iwork)
      !deallocate(ifail)
       CALL DALLOC (ap, 'AP', 'INIEST2')
       CALL DALLOC (eigval, 'EIGVAL', 'INIEST2')
       CALL DALLOC (vec, 'VEC', 'INIEST2')
       CALL DALLOC (work, 'WORK', 'INIEST2')
       CALL DALLOC (iwork, 'IWORK', 'INIEST2')
       CALL DALLOC (ifail, 'IFAIL', 'INIEST2')

      RETURN
      END SUBROUTINE INIEST2
