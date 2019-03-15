!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION GCO (K, IR, IA, IB)
!                                                                      *
!   This routine evaluates a coefficient                               *
!                                                                      *
!                                K                                     *
!                               g   (IA,IB)                            *
!                                IR                                    *
!                                                                      *
!                                                                      *
!   Here  K  is the multipolarity, IR  is the sequence number of the   *
!   configuration, and  IA and IB are orbital sequence  numbers. See   *
!   I P Grant,  B J McKenzie,  P H Norrington,  D F  Mayers, and N C   *
!   Pyper, Computer Phys Commun 21 (1980) 207-231, Eq (7).             *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
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
      INTEGER :: IQA, IQB
      REAL(DOUBLE) :: FAC
!      LOGICAL :: FULLA, FULLB
!-----------------------------------------------
!
         IQA = IQ (IA,IR)
         IQB = IQ (IB,IR)
!GG      IQA = IBITS(IIQA((IA - 1)/4 + 1,IR),8*MOD(IA - 1,4),8)
!GG      IQB = IBITS(IIQA((IB - 1)/4 + 1,IR),8*MOD(IB - 1,4),8)
      IF (IQA==NKJ(IA) + 1 .OR. IQB==NKJ(IB)+1) THEN
         FAC = CLRX(NAK(IA),K,NAK(IB))
         GCO = -DBLE(IQA*IQB)*FAC*FAC
      ELSE
         GCO = 0.0D00
      ENDIF

      RETURN
      END FUNCTION GCO
