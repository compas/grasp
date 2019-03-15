!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION FCO (K, IR, IA, IB)
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               f   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA  and  IB  are orbital  sequence  numbers.   *
!   ( I P Grant,  B J McKenzie,  P H Norrington, D F Mayers, and N C   *
!   Pyper,  Computer Phys Commun 21 (1980) 207-231, Eqs (6). )         *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX, IQ.                                     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
!                                                                      *
!***********************************************************************
! . Two calls to IQ are physically inlined here to reduce the overhead.
! . In the inlined text, the checks for the integer range (ref.
!    lib92/iq.f for details) are removed since FCO is called only by
!    SETCOF and SETHAM where the input arguments always fall in the
!    suitable range.
!XHH 1997.03.05
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE debug_C
      USE orb_C, IIQA=>IQA
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE clrx_I
      USE IQ_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K
      INTEGER  :: IR
      INTEGER , INTENT(IN) :: IA
      INTEGER , INTENT(IN) :: IB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQA, IQF, KAPPA
      REAL(DOUBLE) :: FAC
!-----------------------------------------------
!
      IF (IA == IB) THEN
!
         IQA = IQ (IA,IR)
!         IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
!
         IF (K == 0) THEN
            FCO = DBLE((IQA*(IQA - 1))/2)
         ELSE
            IQF = NKJ(IA) + 1
            IF (IQA == IQF) THEN
               KAPPA = NAK(IA)
               FAC = CLRX(KAPPA,K,KAPPA)*DBLE(IQA)
               FCO = -0.5D0*FAC*FAC
            ELSE
               FCO = 0.D0
            ENDIF
         ENDIF
!
      ELSE
!
         IF (K == 0) THEN
            FCO = DBLE (IQ (IA,IR)*IQ (IB,IR))
!            FCO = DBLE(IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)*IBITS(&
!               IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8))
         ELSE
            FCO = 0.D0
         ENDIF
!
      ENDIF
      RETURN
!*
!  300 FORMAT (/'  ',1I2
!     :        /' f    (',1I2,1A2,',',1I2,1A2,') = ',1PD21.14,
!     :        /'  ',1I3/)
!*
      END FUNCTION FCO
