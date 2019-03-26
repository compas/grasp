      SUBROUTINE INIEST(N, NB, NIV, HMX, JCOL, IROW, BASIS, IBLOCK, JBLOCK)
!-----------------------------------------------------------------------
!     Routine for providing initial estimates from the diagonal
!     of the matrix. This way was used by Dvdson in atomic structure
!       calculations. It should be used to obtain estimates when nothing
!       else is available.
!
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:15:59   2/21/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dinit_I
      !USE dspevx_I
      !USE vec_I
      !USE dcopy_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: NB
      INTEGER  :: NIV
      INTEGER, INTENT(IN) :: JBLOCK
      INTEGER, INTENT(IN) :: JCOL(0:*)
      INTEGER, INTENT(IN) :: IROW(*)
      INTEGER, INTENT(IN) :: IBLOCK(*)
      REAL(DOUBLE), INTENT(IN) :: HMX(*)
      REAL(DOUBLE)  :: BASIS(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: NS, ICOUNT, &
         J, ISTART, ISH, II, NFOUND, INFO, IERR, IC
      !REAL(DOUBLE) :: EIGVAL, WORK

      integer, dimension(:), pointer :: iqiwork,iqif,iqjsh
      integer, dimension(:), allocatable, target :: iwork, ifail, jsh
      real(double), dimension(:), pointer :: iqap,iqeig,iqvec,iqwork
      real(double), dimension(:), allocatable, target :: ap, eigval, vec, work

!-----------------------------------------------
!     !pointer (iqap,ap(1)),(iqeig,eigval(1)),(iqvec,vec(1))
!     !pointer (iqwork,work(1)),(iqiwork,iwork(1)),(iqif,IFAIL(1))
!     !pointer (iqjsh,jsh(1))

!******************************************************************

!***** alloc space for 100*100 lower triangular
      NS = MIN(1000,NB)
!once here        NS = min(800, NB)
        !print*, 'NS=',NS,' NIV=',NIV
      !CALL ALLOC (IQAP, NS*(NS + 1)/2, 8)
      allocate (ap(NS*(NS + 1)/2),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,NS*(NS + 1)/2,'iniest::ap',ierr);
      iqap=>ap

      CALL DINIT (NS*(NS + 1)/2, 0.D0, AP, 1)

      !CALL ALLOC (IQEIG, NS, 8)
      allocate (eigval(NS),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::eigval',ierr);
      IQEIG=>eigval

      !CALL ALLOC (IQVEC, NS*NIV, 8)
      allocate (vec(NS*NIV),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns*niv,'iniest::vec',ierr);
      IQVEC=>vec

      !CALL ALLOC (IQWORK, 8*NS, 8)
      allocate (work(8*NS),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,8*NS,'iniest::work',ierr);
      IQWORK=>work

      !CALL ALLOC (IQIWORK, 5*NS, 4)
      allocate (iwork(5*NS),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,5*NS,'iniest::iwork',ierr);
      IQIWORK=>iwork

      !CALL ALLOC (IQIF, NS, 4)
      allocate (ifail(ns),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::IFAIL',ierr);
      iqif=>ifail

      !CALL ALLOC (IQJSH, NS, 4)
      allocate (jsh(NS),stat=ierr)
      if (ierr.ne.0) call mem_fail(6,ns,'iniest::jsh',ierr);
      iqjsh=>jsh

      ICOUNT = 0

!**** separate upper left block of size NS*NS
      DO J = 1, N
         IF (ICOUNT >= NS) EXIT
         IF (IBLOCK(J) /= JBLOCK) CYCLE
         ICOUNT = ICOUNT + 1
         JSH(ICOUNT) = J
         ISTART = JCOL(J-1) + 1
         ISH = J - ICOUNT
  100    CONTINUE
         AP(IROW(ISTART)-ISH+(ICOUNT-1)*(2*NS-ICOUNT)/2) = HMX(ISTART)
!           print*, 'ic = ',icount, ' i=', irow(istart)-ish
         ISTART = ISTART + 1
! check the block structure for zero elements
         IF (ISTART > JCOL(J)) GO TO 102
         ISH = ISH + COUNT(IBLOCK(IROW(ISTART-1):IROW(ISTART))/=JBLOCK)
         IF (ISTART > JCOL(J)) GO TO 102
         IF (IROW(ISTART) - ISH <= NS) GO TO 100
!           goto 100
  102    CONTINUE
      END DO
      CALL DSPEVX ('Vectors also', 'In a range', 'Lower triangular', NS, AP, &
         -1., -1., 1, NIV, 0.D0, NFOUND, EIGVAL, VEC, NS, WORK, IWORK, IFAIL, &
         INFO)
      IERR = -ABS(INFO)
      IF (IERR /= 0) WRITE (6, *) 'iniest ierr =', IERR
!           print '(D14.7,X,I2)', (eigval(i),i, i=1,NIV)

!******************************************************************


!
!       ..Build the Basis.
!
      CALL DINIT (N*NIV, 0.D0, BASIS, 1)
      ISTART = 0
      DO J = 1, NIV
!scatter the vectors
         DO IC = 1, NS
            BASIS(ISTART+JSH(IC)) = VEC((J - 1)*NS + IC)
         END DO
!              call dcopy(NS,vec((J-1)*NS+1),1,BASIS(istart+ish(j)),1)
         ISTART = ISTART + N
      END DO
      CALL DCOPY (NIV, EIGVAL, 1, BASIS(NIV*N+1), 1)
      deallocate(ap) !CALL DALLOC (IQAP)
      deallocate(eigval) !CALL DALLOC (IQEIG)
      deallocate(vec) !CALL DALLOC (IQVEC)
      deallocate(work) !CALL DALLOC (IQWORK)
      deallocate(iwork) !CALL DALLOC (IQIWORK)
      deallocate(ifail) !CALL DALLOC (IQIF)
      deallocate(jsh) !CALL DALLOC (IQJSH)

      RETURN
      END SUBROUTINE INIEST
