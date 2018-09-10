!***********************************************************************
!
      SUBROUTINE INIESTSD (nmax, ncf, myid, nprocs,                    &
          NIV, BASIS, IMCDF, EAV)
!
!    Routine for providing initial estimates from upper left corner
!       of the matrix. It is exact (not estimates) if ncf <= nmax which
!    is currently set to be 400 in the calling routine.
!
!    Matrix is sparse and on the disk
!
!   Block version by Xinghong He            Last revision: 14 Dec 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dspevx_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER, INTENT(IN) :: nmax, ncf, myid, nprocs, niv, imcdf
       REAL(DOUBLE), INTENT(IN) :: EAV
       REAL(DOUBLE), DIMENSION(*), INTENT(IN):: Basis
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       REAL(DOUBLE), DIMENSION(:), pointer :: ap, eigval, vec, work, hmx
       INTEGER, DIMENSION(:), pointer :: iwork, ifail, irow
       INTEGER :: in, ns, ncfdum, iccutdum, myiddum, nprocsdum, info, &
            ierr, j, joff, nelc, ir, nfound
       REAL(DOUBLE) :: elsto
!-----------------------------------------------------------------------
      NS = min (nmax, ncf)
 
      CALL alloc (ap, (NS*(NS+1))/2, 'AP', 'INIESTSD')
      CALL dinit ((NS*(NS+1))/2, 0.d0, ap, 1)
 
!**** separate upper left block of size NS*NS
 
      CALL alloc (hmx, ncf, 'HMX', 'INIESTSD')
      CALL alloc (irow, ncf, 'IROW', 'INIESTSD')
      READ (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
      IF (myid .EQ. 0) PRINT *, 'iniestsd ...........'
      IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum                       &
                    .OR. nprocsdum .NE. nprocs)                         &
         STOP 'iniestsd: ncf read wrong'
 
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
 
      DO j = myid + 1, ns, nprocs
         joff = (j*(j-1))/2
         READ (IMCDF) NELC,ELSTO,(HMX(IR),IR=1,NELC),                   &
                                (IROW(IR),IR=1,NELC)
         HMX(NELC) = HMX(NELC) - EAV ! Shift the diagonal
         DO ir = 1, nelc
            ap(irow(ir) + joff) = hmx(ir)
         ENDDO
      ENDDO
 
! Let each node have a complete copy of ap
 
!      CALL gdsummpi (ap, (NS*(NS+1))/2)
 
! To be in step with other cases, go through the whole block.
!
! This is not necessary since currently the file pointer is moved
! to the absolute position of the .res files which is always counted
! from the begining of the .res files of each node. Besides, the
! following segment seems not working properly for the last block.
! Xinghong He 98-12-14
 
      !mylast = j - nprocs
      !DO j = mylast, ncf, nprocs
      !   READ (imcdf)
      !ENDDO
 
      CALL dalloc (hmx,'HMX', 'INIESTSD')
      CALL dalloc (irow, 'HMX', 'INIESTSD')
 
      CALL alloc (eigval,NS,'EIGVAL','INIESTSD')
      CALL alloc (vec,NS*NIV,'VEC','INIESTSD')
      CALL alloc (work,8*NS,'WORK','INIESTSD')
      CALL alloc (iwork,5*NS,'IWORK','INIESTSD' )
      CALL alloc (ifail,NS, 'IFAIL','INIESTSD' )
 
      CALL DSPEVX ('Vectors also','In a range','Upper triangular',     &
                NS,AP,-1.0D0,-1.0D0,1,NIV,0.d0,                        &
                NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      IERR = -ABS (INFO)
 
!******************************************************************
 
!       ..Build the Basis.
 
      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)
!       ...scatter the vectors
      DO J = 1, NIV
         CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO
      CALL dcopy (NIV, EIGVAL,1,BASIS(NIV*ncf+1),1)
 
      CALL dalloc (ap, 'AP', 'INIESTSD')
      CALL dalloc (eigval,'EIGVAL', 'INIESTSD')
      CALL dalloc (vec, 'VEC', 'INIESTSD')
      CALL dalloc (work, 'WORK', 'INIESTSD')
      CALL dalloc (iwork, 'IWORK', 'INIESTSD')
      CALL dalloc (ifail, 'IFAIL', 'INIESTSD')
 
      RETURN
      END SUBROUTINE INIESTSD
