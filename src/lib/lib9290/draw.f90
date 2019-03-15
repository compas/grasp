!***********************************************************************
!                                                                      *
      SUBROUTINE DRAW(P, SP, Q, SQ, MF)
!                                                                      *
!   This  subroutine  generates  a  printer plot. P and Q are radial   *
!   functions with the maximum tabulation point MF. SP is the factor   *
!   by which P is to be scaled, SQ is the factor by which Q is to be   *
!   scaled.                                                            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:22   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE GRID_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: MF
      REAL(DOUBLE), INTENT(IN) :: SP
      REAL(DOUBLE), INTENT(IN) :: SQ
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: P
      REAL(DOUBLE), DIMENSION(NNNP), INTENT(IN) :: Q
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IOFFST, IXLOCQ, IXLOCP, IYLOC
      REAL(DOUBLE) :: SIXTY, FIFT4, OHTHT, DMX, DMN, SPI, SQI, RMX, XSCAL, &
         YSCAL, SQX, SPX
      LOGICAL :: FIRST
      CHARACTER :: CBLANK*132, CDASH*132
!cjb
!     CHARACTER , DIMENSION(60) :: CPLOT*132
      CHARACTER , DIMENSION(MF) :: CPLOT*132
!cjb
      CHARACTER :: CP, CQ
!-----------------------------------------------
      DATA SIXTY/ 60.0D00/
      DATA FIFT4/ 54.0D00/
      DATA OHTHT/ 131.0D00/
!
      DATA FIRST/ .TRUE./


!
!   This initialization is carried out once per run
!
      IF (FIRST) THEN
         CBLANK(1:1) = '|'
         CDASH(1:1) = '-'
         DO I = 2, 131
            CBLANK(I:I) = ' '
            CDASH(I:I) = '-'
         END DO
         CBLANK(132:132) = '|'
         CDASH(132:132) = '-'
         FIRST = .FALSE.
      ENDIF
!
!   Initialization
!
      IF (SQ == 0.0D00) THEN
         CP = 'X'
         CQ = '-'
      ELSE
         CP = 'P'
         CQ = 'Q'
      ENDIF
!
!   Determine the range of amplitude
!
      DMX = MAX(SP*P(1),SQ*Q(1))
      DMN = MIN(SP*P(1),SQ*Q(1))
      DO I = 2, MF
         SPI = SP*P(I)
         SQI = SQ*Q(I)
         DMX = MAX(MAX(DMX,SPI),SQI)
         DMN = MIN(MIN(DMN,SPI),SQI)
      END DO
!
!   Determine the radial extent of the function
!
      RMX = R(MF)
!
!   Determine the scale factors
!
      IF (DMX==0.0D00 .AND. DMN==0.0D00) THEN
         WRITE (99, 300)
         RETURN
      ELSE
         XSCAL = FIFT4/ABS(DMX - DMN)
      ENDIF
      YSCAL = OHTHT/DBLE(N)
!
!   Locate x = 0 if this is in the range
!
      IF (DMX>0.0D00 .AND. DMN<0.0D00) THEN
         IOFFST = 4 + XSCAL*ABS(DMN)
      ELSE
         IOFFST = 1
      ENDIF
!
!   Initialize the array CPLOT
!
      DO I = 1, 60
         IF (I==1 .OR. I==IOFFST .OR. I==60) THEN
            CPLOT(I) = CDASH
         ELSE
            CPLOT(I) = CBLANK
         ENDIF
      END DO
!
!   Generate the plot
!
!   Note that 'P' is the dominant character
!
      SQX = SQ*XSCAL
      SPX = SP*XSCAL
      DO I = 1, MF
         IXLOCQ = NINT(SQX*Q(I)) + IOFFST
         IXLOCP = NINT(SPX*P(I)) + IOFFST
         IYLOC = MAX(1,NINT(DBLE(I)*YSCAL))
         CPLOT(IXLOCQ)(IYLOC:IYLOC) = CQ
         CPLOT(IXLOCP)(IYLOC:IYLOC) = CP
      END DO
!
!   Print plot
!
      WRITE (99, 301)
      DO I = 60, 1, -1
         WRITE (99, 302) CPLOT(I)
      END DO
!
!   Plot-information line
!
      WRITE (99, 303) RMX, DMX
!
      RETURN
!
  300 FORMAT(' No plot: function is identically zero')
  301 FORMAT('1')
  302 FORMAT(1X,1A132)
  303 FORMAT(/,50X,1P,' r (max) = ',1D10.3,'Bohr radii,',&
         ' Maximum of functions = ',1D10.3)
      RETURN
!
      END SUBROUTINE DRAW
