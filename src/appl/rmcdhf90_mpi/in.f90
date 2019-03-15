!***********************************************************************
!                                                                      *
      SUBROUTINE IN(IORB, JP, P, Q, MTP)
!                                                                      *
!   This program computes the solution of an  inhomogeneous  pair of   *
!   radial Dirac equations in the tail region. A simple extension of   *
!   the method of  C Froese Fischer [C Froese, Can J Phys, 41 (1963)   *
!   1895]  is  used. The  equations  are treated as a boundary value   *
!   problem, with  the value of P(r) given at the inner boundary and   *
!   required to  be  sufficiently small  at the outer boundary.  The   *
!   location  of this second boundary is determined in the course of   *
!   the calculation. The same  finite  difference equations are used   *
!   as in the outward integration.                                     *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      IORB: (Input) Index of orbital                                  *
!      JP:   (Input) Tabulation point at which the outward             *
!            integration terminated                                    *
!      P,Q:  (Input and Output) arrays containing, respectively, the   *
!            large and small components of the solution                *
!      MTP:  (Output) Maximum tabulation point of functions P, Q       *
!                                                                      *
!   When  written  in  matrix form the system of linear equations iS   *
!   M*W = V, where  W  is a vector consisting of the elements Q(JP),   *
!   P(JP+1), Q(JP+1), ... . The matrix  M  is  of band type, with  5   *
!   elements  in each row.  M is  expressed as a product of two tri-   *
!   angular matrices L and  U (the L/U decomposition, carried out by   *
!   the Crout procedure), variable elements of L being stored in the   *
!   arrays TC, TD  and TE, and  variable elements of U in TH, TI and   *
!   TJ. The  solution vector W  is then obtained by solving two tri-   *
!   angular systems, L*Z = V AND U*W = Z.                              *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                           Last update: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNP
      USE def_C,            ONLY: accy
      USE grid_C,           ONLY: h, n, r, rpor
      USE int_C,            ONLY: TF, TG, XU, XV
      USE orb_C,            ONLY: nak
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE convrt_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IORB
      INTEGER , INTENT(IN) :: JP
      INTEGER , INTENT(OUT) :: MTP
      REAL(DOUBLE) , INTENT(INOUT) :: P(NNNP)
      REAL(DOUBLE) , INTENT(INOUT) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, LCNUM, MTPP1
      REAL(DOUBLE), DIMENSION(NNNP) :: TH, TI, TJ, XR, XS
      REAL(DOUBLE) :: EPS, HHK, HHKF, CPJ, CMJ, CPJP1, CMJP1, TEI, TTHIS, TCI, &
         TDI, TLAST
      CHARACTER :: CNUM*5
!-----------------------------------------------
!
!   Global initializations
!
!ww      EPS = 0.1D 00*ACCY
      EPS = 0.01D00*ACCY
      HHK = 0.5D00*H*DBLE(NAK(IORB))
!
!   Initialize counters
!
      I = 1
      J = JP
!
!   Other initializations
!
      HHKF = HHK*RPOR(J)
      CPJ = 1.0D00 + HHKF
      CMJ = 1.0D00 - HHKF
      HHKF = HHK*RPOR(J+1)
      CPJP1 = 1.0D00 + HHKF
      CMJP1 = 1.0D00 - HHKF
!
!   Compute required elements of first two rows of  L  and  U
!
      TH(I) = CPJ
      TEI = -TF(J)/TH(I)
      TI(I) = (-CPJP1) + TEI*TG(J+1)
      TJ(I) = CMJP1*TEI - TF(J+1)
!
!   First elements of solution vector Z
!
      XR(I) = (-XV(J)) + TG(J)*P(J)
      XS(I) = (-XU(J)) - CMJ*P(J) - TEI*XR(I)
!
      TTHIS = ABS(XS(I)/TI(I))
!
    1 CONTINUE
      I = I + 1
      J = J + 1
!
!   Failure if tables not long enough
!
      IF (J >= N) THEN
         CALL CONVRT (N, CNUM, LCNUM)
         WRITE (6, *) 'IN: maximum tabulation point exceeds'
         WRITE (6, *) ' dimensional limit (currently '//CNUM(1:LCNUM)//');'
         WRITE (6, *) ' radial wavefunction may indicate a'
         WRITE (6, *) ' continuum state.'
         STOP
      ENDIF
!
!   Compute required elements of remaining rows of  L  and  U
!
      CPJ = CPJP1
      CMJ = CMJP1
      HHKF = HHK*RPOR(J+1)
      CPJP1 = 1.0D00 + HHKF
      CMJP1 = 1.0D00 - HHKF
!
      TCI = -TG(J)/TI(I-1)
      TH(I) = CPJ - TCI*TJ(I-1)
      TDI = CMJ/TI(I-1)
      TEI = ((-TF(J))-TDI*TJ(I-1))/TH(I)
      TI(I) = (-CPJP1) + TEI*TG(J+1)
      TJ(I) = (-TF(J+1)) + CMJP1*TEI
!
!   Solution of  L*Z = V
!
      XR(I) = (-XV(J)) - TCI*XS(I-1)
      XS(I) = (-XU(J)) - TDI*XS(I-1) - TEI*XR(I)
!
!   Test for outer boundary
!
      TLAST = TTHIS
      TTHIS = ABS(XS(I)/TI(I))
      IF (TTHIS + TLAST <= EPS) THEN
         MTP = J
      ELSE
         GO TO 1
      ENDIF
!
!   Reset counter
!
      I = I - 1
!
!   Last two rows of solution of  U*W = Z ; evaluation of Q(J)
!
      Q(J) = 0.0D00
      P(J) = XS(I)/TI(I)
      Q(J-1) = (XR(I)+TG(J)*P(J))/TH(I)
!
!   Solution of  U*W = Z
!
    2 CONTINUE
      J = J - 1
      I = I - 1
!
      IF (I > 0) THEN
         P(J) = (XS(I)-TJ(I)*Q(J))/TI(I)
         Q(J-1) = (XR(I)+(1.0D00-HHK*RPOR(J))*Q(J)+TG(J)*P(J))/TH(I)
         GO TO 2
      ENDIF
!
!   Complete tables with zeroes
!
      MTPP1 = MTP + 1
      P(MTPP1:N) = 0.0D00
      Q(MTPP1:N) = 0.0D00
!
      RETURN
      END SUBROUTINE IN
