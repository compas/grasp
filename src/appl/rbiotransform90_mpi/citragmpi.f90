!***********************************************************************
!                                                                      *
!                                                                      *
      SUBROUTINE CITRAG(CIIN,NCSF,NCIV,L,NSHL,T,NIN,NTESTG,CIOUT,SCR)
!                                                                      *
!   Calculate the action of the  operator                              *
!                                                                      *
!    T(L) = T(nl) T(n-1 L ) .... T(1 l)                                *
!                                                                      *
!   where                                                              *
!                                                                      *
!    T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))     *
!    t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)                             *
!    E(lj,li) total symmetric shell excitation operator                *
!                                                                      *
!   Modification of the CITRA routine for MCHF                         *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN      : Input CI vectors ( Destroyed in the process )          *
!   NCSCF     : Length of CI expansion                                 *
!   NCIV      : Number of CI vectors                                   *
!   L         : L value of shells in excitations                       *
!   NSHL      : Number of shells with this L                           *
!   T         : List of coefficients                                   *
!   NIN       : Number of inactive shells                              *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE sbdat1_C
      USE orb_C
      USE mpi_C
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE wrtmat_I
      USE scalve_I
!      USE copvec_I
      USE tiinig_I
      USE ti1tv_I
      USE vecsum_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCSF
      INTEGER  :: NCIV
      INTEGER  :: L
      INTEGER  :: NSHL
      INTEGER, INTENT(IN) :: NIN
      INTEGER  :: NTESTG
      REAL(DOUBLE)  :: CIIN(NCSF,NCIV)
      REAL(DOUBLE)  :: T(NSHL,NSHL)
      REAL(DOUBLE)  :: CIOUT(NCSF,NCIV)
      REAL(DOUBLE)  :: SCR(NCSF,NCIV)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, IIN, IPOT, I, N
      REAL(DOUBLE) :: FACTOR, TII, XNFACI
      REAL(DOUBLE)  :: SCRtmp(NCSF,NCIV)
!-----------------------------------------------
!

      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 0

!
      IF (NTEST >= 10) THEN
         WRITE (6, *)
         WRITE (6, *) ' ***************'
         WRITE (6, *) ' Entering CITRAG'
         WRITE (6, *) ' ***************'
         WRITE (6, *)
      ENDIF
      IF (NTEST >= 100) THEN
         WRITE (6, *) ' Input CI vectors '
         CALL WRTMAT (CIIN, NCSF, NCIV, NCSF, NCIV)
         WRITE (6, *) ' Transformation matrix T'
         CALL WRTMAT (T, NSHL, NSHL, NSHL, NSHL)
      ENDIF
!
!. Factor from inactive shells
!
      IF (NIN /= 0) THEN
         FACTOR = 1.0D0
         DO IIN = 1, NIN
            FACTOR = FACTOR*T(IIN,IIN)
         END DO
!
!       IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
!                          by (2j+1) corresponding to L)
!
         IPOT = 2*IABS(NAK(NSHLP(L,IIN)))
         FACTOR = FACTOR**IPOT
         CALL SCALVE (CIIN, FACTOR, NCIV*NCSF)
      ENDIF
      IF (NIN == NSHL) CALL COPVEC (CIIN, CIOUT, NCIV*NCSF)
!
      DO I = NIN + 1, NSHL
         IF (NTEST >= 100) WRITE (6, *) ' Loop I,L = ', I, L
!
!. The diagonal contribution
!
         TII = T(I,I)
         CALL TIINIG (CIIN, NCSF, NCIV, I, L, TII, CIOUT, NTESTG)
!        IF (LWORK2.LT.NCIV*NCSF) THEN
!          WRITE(*,*) 'In CITRAG: Dimension of LWORK2 must be',
!     &               'increased to at least',NCIV*NCSF
!        ENDIF
         CALL COPVEC (CIOUT, SCR, NCIV*NCSF)
!
!.  Off diagonal contributions
!
         XNFACI = 1.0D0
         DO N = 1, 2*IABS(NAK(NSHLP(L,I)))
            IF (NTEST >= 100) WRITE (6, *) ' Loop N = ', N
!
!   T ** (N-1) is supposed to be in SCR, copy to CIIN
!   and apply S
!
            CALL COPVEC (SCR, CIIN, NCIV*NCSF)
            CALL TI1TV (CIIN,NCSF,NCIV,I,L,T(1,I),NSHL,SCRtmp,NTESTG)
!            CALL MPI_ALLREDUCE(SCRtmp(1,1),SCR(1,1),NCIV*NCSF,         &
            CALL MPI_ALLREDUCE(SCRtmp,SCR,NCIV*NCSF,         &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            XNFACI = XNFACI/FLOAT(N)
            CALL VECSUM (CIOUT, CIOUT, SCR, 1.0D0, XNFACI, NCIV*NCSF)
         END DO
         CALL COPVEC (CIOUT, CIIN, NCIV*NCSF)
!
      END DO
!
      IF (NTEST >= 100) THEN
         WRITE (6, *) ' Output CI vectors L = ', L
         CALL WRTMAT (CIOUT, NCSF, NCIV, NCSF, NCIV)
      ENDIF
!
      RETURN
      END SUBROUTINE CITRAG
