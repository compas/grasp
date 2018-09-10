!***********************************************************************
!                                                                      *
      SUBROUTINE MANEIG(dvdfirst, LPRINT, JBLOCK, &
                        NCFPAT, NCMINPAT, NEVECPAT, NCFTOT) 
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: eigenpairs are found        *
!          using DVDSON; this involves up to three steps:              *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                                                                      *
!   We  assume that  the sparse representation  of the matrix  is in   *
!   core.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, RALLOC.          *
!               [RSCF92]: POSNFL, SPICMV2.
!               [DVDSON]: DVDSON                                       *
!               [AUXBLAS]: DINIT/SINIT                                 *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified by Xinghong He               Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE memory_man
      USE DEF_C 
      USE eigv_C
      USE hblock_C
      USE hmat_C
      USE mpi_s
      USE orb_C
      USE WCHBLK_C, JBLOCKK=>JBLOCK 
      USE WHERE_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE spicmv2_I 
      USE iniest2_I 
      USE gdvd_I 
      USE itjpo_I 
      USE ispar_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, INTENT(IN) ::  dvdfirst
      INTEGER  :: JBLOCK 
      INTEGER, INTENT(IN) :: NCFPAT 
      INTEGER, INTENT(IN) :: NCMINPAT 
      INTEGER, INTENT(IN) :: NEVECPAT 
      INTEGER  :: NCFTOT 
      LOGICAL  :: LPRINT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NVECT, LIM, LWORK, NVEX, NIV, MAXITR, N1000, MBLOCK, &
         ILOW, IHIGH, LIWORK, IC, NLOOPS, NMV, NEND, J, JSTATE, &
         IOFSET, I, IA, IERR
      INTEGER, DIMENSION(:), POINTER :: IWORK, JWORK
      REAL(DOUBLE) :: PNWORK, CRITE, CRITC, CRITR, ORTHO, AMAX, WA, DNFAC 
      REAL(DOUBLE), DIMENSION(:), POINTER :: WORK, DIAG
      LOGICAL :: HIEND 
!-----------------------------------------------

      !PRINT *, 'maneig ...'
 
!      ...spicmv2 needs this COMMON /WCHBLK/JBLOCKK
      JBLOCKK = JBLOCK 
!
!=======================================================================
!   Trivial case
!=======================================================================
      IF (NCF == 1) THEN 
         EVAL(NCMINPAT+1) = 0.D0 
         EVEC(NEVECPAT+1) = 1.D0 
         GO TO 123                               ! Don't like big ELSE 
      ENDIF 
!
!=======================================================================
!   Non-trivial case - Use Davidson eigensolver
!=======================================================================
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECT
!   plus a smaller number if storage is severely constrained
!
      NVECT = NCMAXBLK(JBLOCK) 
      LIM = MIN(NCF,2*NVECT + 40) 
      LWORK = 2*NCF*LIM + LIM*LIM*2 + 11*LIM + NVECT 
      CALL ALLOC (WORK, LWORK, 'WORK', 'MANEIG') 
 
      !...At most 14 ? restriction removed xhh 98-05-19
      !nvex = MIN (nvect,ncfblk(jblock),14)
      NVEX = MIN(NVECT,NCFBLK(JBLOCK)) 
      NIV = NVEX 
      MAXITR = MIN(NVECT*200,NCF) 
!      N1000 = 2000 
      N1000 = 4000 
!
!   Initial estimates for eigenvectors
!
!CFF
      if (dvdfirst .or. (ncf .LE. n1000) ) then
         CALL INIEST2 (N1000, NCF, NIV, WORK, EMT, IENDC, IROW) 
      else
!CFF        .. use current estimates
      nend = ncf * nvex
         DO j = 1, nevblk(jblock)
            work( nend + iccmin(j+ncminpat) ) = eval(ncminpat+j)
            CALL dcopy (ncf, evec( nevecpat + ncf*(j-1) + 1)    ,    1, &
                          work( ncf*(iccmin(j+ncminpat)-1) + 1 ), 1)
         ENDDO
      ENDIF
 
! iniest2 looks for eigenvectors of n1000*n1000 matrix so there
! is no need to call dvdson if block size <= n1000
 
      IF (NCF > N1000) THEN 
          WRITE (*, *) 'Calling dvdson!!!', MAXITR 
 
         ! Call Davidson eigensolver
 
         MBLOCK = 1 
         ILOW = 1 
         IHIGH = NVEX 
         LIWORK = 6*LIM + NVECT 
         CRITE = 1.0D-17 
!         CRITC = 1.0D-08 
!         CRITR = 1.0D-08 
!         ORTHO = MAX(1D-8,CRITR) 
         critc = 1.0D-09
         critr = 1.0D-09
         ortho = MAX (1D-9, critr)
!
!   Store the diagonals in a separate array and make it global
!
         CALL ALLOC (DIAG, NCF, 'DIAG', 'MANEIG') 
 
         DO IC = MYID + 1, NCF, NPROCS 
            DIAG(IC) = EMT(IENDC(IC)) 
         END DO 
 
         CALL ALLOC (IWORK, LIWORK, 'IWORK', 'MANEIG') 
         CALL ALLOC (JWORK, LIM, 'JWORK', 'MANEIG') 
         CALL GDVD (SPICMV2, NCF, LIM, DIAG, ILOW, IHIGH, JWORK, NIV, MBLOCK, &
            CRITE, CRITC, CRITR, ORTHO, MAXITR, WORK, LWORK, IWORK, LIWORK, &
            HIEND, NLOOPS, NMV, IERR) 
 
         CALL DALLOC (DIAG, 'DIAG', 'MANEIG') 
         CALL DALLOC (IWORK, 'IWORK', 'MANEIG') 
         CALL DALLOC (JWORK, 'JWORK', 'MANEIG') 
 
         WRITE (*, 301) NLOOPS, NMV 
         IF (IERR /= 0) WRITE (*, 302) IERR 
      ENDIF 
!
!   Pick up the eigen pairs and store in EVAL and EVEC
!
      NEND = NCF*NVEX 
      DO J = 1, NEVBLK(JBLOCK) 
         EVAL(NCMINPAT+J) = WORK(NEND + ICCMIN(J + NCMINPAT)) 
         CALL DCOPY (NCF, WORK(NCF*(ICCMIN(J + NCMINPAT) - 1) + 1), 1, EVEC(&
            NEVECPAT+NCF*(J-1)+1), 1) 
      END DO 
!     print *, ncminpat,(eval(ncminpat+j),j=1,nevblk(jblock)),
!    1 'zou,from maneig'
!
!   Deallocate storage
!
      CALL DALLOC (WORK, 'WORK', 'MANEIG') 
 
  123 CONTINUE 
      DO JSTATE = 1, NEVBLK(JBLOCK) 
!
!   Find the dominant component of each eigenvector
!
         IOFSET = NEVECPAT + NCF*(JSTATE - 1) 
 
         AMAX = 0.D0 
         DO I = 1, NCF 
            WA = ABS(EVEC(I+IOFSET)) 
            IF (WA <= AMAX) CYCLE  
            AMAX = WA 
            IA = I 
         END DO 
!
!   Find the angular momentum and parity of the dominant component
!
!GG         IATJPO(JSTATE+NCMINPAT) = ITJPO(IA + NCFPAT) 
!GG         IASPAR(JSTATE+NCMINPAT) = ISPAR(IA + NCFPAT) 
!
!   Redefine eigenvectors so that the dominant component
!   is positive
!
         IF (EVEC(IA+IOFSET) >= 0.D0) CYCLE  
         DNFAC = -1.D0 
         CALL DSCAL (NCF, DNFAC, EVEC(IOFSET+1), 1) 
!===============================================================
 
      END DO 
 
  301 FORMAT('DVDSON: ',1I3,' loops; ',1I3,' matrix-vector multiplies.') 
  302 FORMAT(' Returned from DVDSON with IERR = ',1I4) 
  303 FORMAT(/,' ***** WARNING *****'/,/,&
         ' The angular momentum and parity of level ',1I2,' have changed:'/,&
         ' Last iteration: (2J+1) = ',1I2,', parity = ',1I2,';'/,&
         ' this iteration: (2J+1) = ',1I2,', parity = ',1I2,'.') 
 
      RETURN  
      END SUBROUTINE MANEIG 
