!***********************************************************************
!                                                                      *
      SUBROUTINE WGHTD5(iatjpo, iaspar)
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: ALLOC, DALLOC, ISPAR, ITJPO.                           *
!                                                                      *
!                                          Last updated: 02 Nov 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE eigv_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE orb_C, ONLY: ncf, nw, iqa
      USE prnt_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE ispar_I
!      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: iatjpo, iaspar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER, DIMENSION(5) :: ICONF
!      REAL(DOUBLE), DIMENSION(5) :: WGHT(5)
      REAL(DOUBLE), DIMENSION(5) :: WT(5)
!      REAL(DOUBLE), DIMENSION(:), pointer :: wt
!      INTEGER, DIMENSION(:), pointer :: next
      REAL(DOUBLE) :: awt
      INTEGER :: i, j, jj, ibegin, iend, nelt, iv, icf, ip
!
!   Allocate storage for local arrays
!
!CFF  These arrays are not needed
!      CALL ALLOC (WT, NCF,'WT', 'WGHTD5')
!      CALL ALLOC (NEXT,NCF,'NEXT', 'WGHTD5')
!
      WRITE (24,300)
!
!CFF  Initialized NELT for NCF >5
      NELT = 5
      IF (NCF .LT. 5) NELT = NCF
      DO 7 IV = 1,NVEC
!
         ICF = IVEC(IV)
!        .., initialize for a new eigenvector
         wt = 0.d0;  iconf = 0
         ibegin = (IV-1)*NCF
         iend   = IV*NCF
!
         DO 4 I = ibegin+1, iend
            awt = abs(EVEC(i))
            if (awt > wt(5)) then
!              ... the weight needs to be inserted
               do j = 1,5
                  if ( awt > wt(j) ) then
                     if (j < 5) then
!                       ... shift right
                        do jj = 5, j+1,-1
                           wt(jj) = wt(jj-1)
                           iconf(jj) = iconf(jj-1)
                        end do
                     end if
!                    ...  insert the element
                     wt(j) = awt
                     iconf(j) = i
                     exit
                  end if
               end do
            end if
    4    CONTINUE
!
!   Print first five elements of list.
!
         IP = (iaspar + 3 ) / 2
!        NELT = MIN (I,5)
         WRITE (24,301) ICF,LABJ(IATJPO),LABP(IP),      &
               (EVEC(iconf(j)),j = 1,NELT)
         WRITE (24,302) (ICONF(j)-ibegin,j = 1,NELT)
    7 CONTINUE
!
!   Deallocate storage for local arrays
!
!      CALL DALLOC (WT, 'WT', 'WGHTD5')
!      CALL DALLOC (NEXT, 'NEXT', 'WGHTD5')
!
      RETURN
!
  300 FORMAT (/'Weights of major contributors to ASF:'                  &
             //'Level J Parity      CSF contributions'/)
  301 FORMAT (I3,2X,2A4,5(3X,F8.5))
  302 FORMAT (13X,5(I11))
!
      END SUBROUTINE WGHTD5
