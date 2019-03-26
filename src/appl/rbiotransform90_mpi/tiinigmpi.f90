!***********************************************************************
!                                                                      *
      SUBROUTINE TIINIG(CIIN, NCSF, NCIV, I, L, CONST, CIOUT, NTESTG)
!                                                                      *
!   Calculates the action of the operator                              *
!   Const ** E(li,li) on a set of vectors                              *
!                                                                      *
!   Adapted for GRASP, daughter of TIINI, born February 1996           *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN      : Input CI vectors                                       *
!   NCSCF     : Length of CI expansion                                 *
!   NCIV      : Number of CI vectors                                   *
!   I         : Shell number                                           *
!   L         : symmetry                                               *
!   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
!               the K:th shell with symmetry L                         *
!   CONST     : The constant                                           *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:41:42   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE mcpdata_C
      USE sbdat1_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setvec_I
      USE wrtmat_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCSF
      INTEGER  :: NCIV
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: NTESTG
      REAL(DOUBLE), INTENT(IN) :: CONST
      REAL(DOUBLE), DIMENSION(NCSF,NCIV), INTENT(IN)    :: CIIN
      REAL(DOUBLE), DIMENSION(NCSF,NCIV), INTENT(INOUT) :: CIOUT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, NFOUND, K, IA, IB, IFIRST, IELMNT, IVAL, ILEFT&
         , IVEC
      REAL(DOUBLE), DIMENSION(NCSF,NCIV) :: CIOUTtmp
      REAL(DOUBLE) :: CONSTN
!-----------------------------------------------
!
!
      NTESTL = 0
      NTEST = MAX(NTESTL,NTESTG)

      IF (NTEST >= 10) WRITE (6, *) ' Entering TIINI'

      CALL SETVEC (CIOUT, 0.0D0, NCSF*NCIV)
      if (NCOEFF .EQ. 0) goto 200
!
!.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
!.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
!.
!     NFOUND : Number of coefficients obtained
!     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
!     ILEFT  = CSF(L) ?
!
      NFOUND = 0
      DO K = 1, NINTG
         IA = INTGRL(K)/KEY
         IB = MOD(INTGRL(K),KEY)
         IF (NSHLP(L,I)/=IA .OR. IA/=IB) CYCLE
         IF (K == 1) THEN
            IFIRST = 1
         ELSE
            IFIRST = INTPTR(K - 1) + 1
         ENDIF
         NFOUND = INTPTR(K) - IFIRST + 1
         EXIT
      END DO

      DO IELMNT = 1, NFOUND
! Bug 2011-08-18 Per Jonsson        IVAL = CNN(IFIRST-1+IELMNT)
!GG         IVAL = IDNINT(CNN(IFIRST-1+IELMNT))
         IVAL = CNN(IFIRST-1+IELMNT)
         CONSTN = CONST**IVAL
         ILEFT = JANN(IFIRST - 1 + IELMNT)
         CIOUT(ILEFT,:NCIV) = CONSTN*CIIN(ILEFT,:NCIV)
      END DO
 200  CONTINUE
      call copvec(CIOUT, CIOUTtmp, NCIV*NCSF)
!GG      call MPI_ALLREDUCE(CIOUTtmp(1,1), CIOUT(1,1), NCIV*NCSF,      &
      call MPI_ALLREDUCE(CIOUTtmp, CIOUT, NCIV*NCSF,      &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!
!.  The previous provided us with all
!   terms with nonvanishing occupation.
!   For terms with vanishing occupation of il,
!   just copy coefficients, since (x) ** 0 = 1
!
      WHERE (CIOUT(:NCSF,:NCIV) == 0.0D0)
         CIOUT(:NCSF,:NCIV) = CIIN(:NCSF,:NCIV)
      END WHERE
!
      IF (NTEST >= 100) THEN
         WRITE (6, *)
         WRITE (6, *) ' Input and output vectors from TIINI I,L', I, L
         CALL WRTMAT (CIIN, NCSF, NCIV, NCSF, NCIV)
         WRITE (6, *)
         CALL WRTMAT (CIOUT, NCSF, NCIV, NCSF, NCIV)
         WRITE (6, *)
      ENDIF
!
      IF (NTEST >= 10) WRITE (6, *) ' Leaving  TIINI'

      RETURN
      END SUBROUTINE TIINIG
