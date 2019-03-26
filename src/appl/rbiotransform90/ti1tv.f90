!***********************************************************************
!                                                                      *
      SUBROUTINE TI1TV(CIIN,NCSF,NCIV,I,L,T,NSHL,CIOUT,NTESTG)
!                                                                      *
!   Calculate the action of the  operator                              *
!   Sum( j) T(jl) E(jl,il) on a set of vectors                         *
!                                                                      *
!   only j<i terms ! included
!                                                                      *
!   Modification of the TI1TV program for MCHF                         *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN  : Input CI vectors                                           *
!   NCSCF : Length of CI expansion                                     *
!   NCIV  : Number of CI vectors                                       *
!   I     : Shell number                                               *
!   L     : L value of shells of this excitation                       *
!   NSHL  : Number of shells with this L                               *
!   T     : The integrals T(*,il)                                      *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY: KEYORB
      USE mcpdata_C
      USE sbdat1_C
      USE orbord_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setvec_I
      USE wrtmat_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ncsf, nciv, nshl, i, L, ntestg
      REAL(DOUBLE), DIMENSION(ncsf, nciv), INTENT(IN) :: ciin
      REAL(DOUBLE), DIMENSION(nshl), INTENT(IN) :: t
      REAL(DOUBLE), DIMENSION(ncsf, nciv), INTENT(OUT) :: ciout
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY=KEYORB
      INTEGER :: ntestL, ntest, nfound, k, ia, ib, ir, ira, irb, ifirst, &
                 j, ileft, iright, ivec, ielmnt
      REAL(DOUBLE) :: racah, x
!
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      NTEST = 000
!
      IF(NTEST.GE.10) WRITE(6,*) ' Entering TI1TV'

      CALL SETVEC(CIOUT,0.0D0,NCIV*NCSF)
!
!.  Obtain address of first coupling coefficient for h(jl,il) :  IFIRST
!.  Obtain number of        coupling coefficient for h(jl,il) :  NFOUND
!
!        NFOUND : Number of coefficients obtained
!
      NFOUND = 0
      DO K = 1,NINTG
        IA = INTGRL(K)/KEY
        IB = MOD(INTGRL(K),KEY)
!
!   Determine the ordernumber of the orbitals I, IA and IB within the
!   shell with symmetry L
!
        IR = NSHLPP(L,NSHLP(L,I))
        IRA = NSHLPP(L,IA)
        IRB = NSHLPP(L,IB)

        IF (IR.EQ.IRB.AND.IRA.LT.IRB) THEN
          IF (K.EQ.1) THEN
            IFIRST = 1
          ELSE
            IFIRST = INTPTR(K-1) + 1
          ENDIF
          NFOUND = INTPTR(K) - IFIRST + 1

   15     DO IELMNT = 1, NFOUND
            RACAH = CNN(IFIRST-1+IELMNT)
            J = IRA
            X = T(J)*RACAH
!ww            IF (NORDII.EQ.0.AND.NORDFF.EQ.0) THEN
!
!  Expression for normal orbital ordering
!
              ILEFT = JANN(IFIRST-1+IELMNT)
              IRIGHT = JBNN(IFIRST-1+IELMNT)
!ww            ELSEIF (NORDII.EQ.1.AND.NORDFF.EQ.1) THEN
!
!  Expression for reverse orbital ordering
!
!ww              IRIGHT = JANN(IFIRST-1+IELMNT)
!ww              ILEFT = JBNN(IFIRST-1+IELMNT)
!ww            ELSE
!ww              WRITE(*,*) 'SOMETHING WRONG'
!ww              STOP
!ww            ENDIF
            DO IVEC =1, NCIV
              CIOUT(ILEFT,IVEC) =                                       &
     &        CIOUT(ILEFT,IVEC) + X*CIIN(IRIGHT,IVEC)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input and output vectors from TI1TV I,L',I,L
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
      END IF
!
      IF(NTEST.GE.10) WRITE(6,*) ' LEAVING TI1TV'

      RETURN
      END
