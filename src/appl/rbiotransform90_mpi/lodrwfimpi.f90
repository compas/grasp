!***********************************************************************
!                                                                      *
      SUBROUTINE LODRWFI(NAME, NTESTG)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:29:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW, NNNP
      USE memory_man
      USE biorb_C
      USE def_C,           ONLY: z, c
      USE DEBUG_C
      USE grid_C
      USE npar_C
      USE sbdat_C,         ONLY: kamax, nshlii, nshlpii
      USE wave_C,          ONLY: pfii, qfii, pzii, mfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE intrpqi_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NTESTG
      CHARACTER, INTENT(IN) :: NAME*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(NNNW) :: NAK, NP
      INTEGER :: NTESTL, NTEST, J, K, I, NWIN, IOS, NPY, NAKY, MY, JJ, KK, JJJ&
         , KKK, IERR
      REAL(DOUBLE), DIMENSION(NNNW) :: E, GAMA
      REAL(DOUBLE) :: CON, FKK, EY, PZY, DNORM
      REAL(DOUBLE), DIMENSION(:), pointer :: pa, qa, ra
      CHARACTER , DIMENSION(NNNW) :: NH*2
      CHARACTER :: G92RWF*6
!-----------------------------------------------
!
!  Common relevant for the initial state
!
!
      NTESTL = 0
      NTEST = MAX0(NTESTL,NTESTG)
      NTEST = 0

!
!   Write entry message
!
      WRITE (6, *) 'Loading Radial WaveFunction File for initial state...'
!
!   Open the radial wave function file
!
      J = INDEX(NAME,' ')
      OPEN(UNIT=21, FILE=NAME(1:J-1)//'.w', FORM='UNFORMATTED', STATUS='OLD', &
         POSITION='asis')
!
!   Save NAK, NP and NH
!
      NAK(:NWII) = NAKII(:NWII)
      NP(:NWII) = NPII(:NWII)
      NH(:NWII) = NHII(:NWII)
!
!   Allocate storage to orbital arrays
!
      CALL ALLOC (PFII, NNNP,NWII, 'PFII', 'LODRWII')
      CALL ALLOC (QFII, NNNP,NWII, 'QFII', 'LODRWII')

!
      CON = Z/C
      CON = CON*CON
!
      DO J = 1, NWII
         PFII(:NNNP,J) = 0.0D00
         QFII(:NNNP,J) = 0.0D00
!
         K = ABS(NAK(J))
         IF (NPARM /= 0) CYCLE
         FKK = DBLE(K*K)
         IF (FKK >= CON) THEN
            GAMA(J) = SQRT(FKK - CON)
         ELSE
            WRITE (6, *) 'LODRWF: Imaginary gamma parameter'
            WRITE (6, *) ' for ', NP(J), NH(J), ' orbital; the'
            WRITE (6, *) ' point model for the nucleus'
            WRITE (6, *) ' is inappropriate for Z > ', C, '.'
            STOP
         ENDIF
!
      END DO
!
!   Read orbital information from Read Orbitals File;
!
      NWIN = 0
      READ (21, IOSTAT=IOS) G92RWF
      IF (IOS/=0 .OR. G92RWF/='G92RWF') THEN
         WRITE (6, *) 'This is not a Radial WaveFunction File;'
         CLOSE(21)
      ENDIF

      IF (NTEST >= 100) THEN
         WRITE (*, *) '******************'
         WRITE (*, *) ' Entering lodrwfi'
         WRITE (*, *) '******************'
      ENDIF

    3 CONTINUE
      READ (21, IOSTAT=IOS) NPY, NAKY, EY, MY
      IF (IOS == 0) THEN
         CALL ALLOC (PA, MY, 'PA', 'LODRWFII')
         CALL ALLOC (QA, MY, 'QA', 'LODRWFII')
         CALL ALLOC (RA, MY, 'RA', 'LODRWFII')
         READ (21) PZY, (PA(I),I=1,MY), (QA(I),I=1,MY)
         READ (21) (RA(I),I=1,MY)
!
!    Orbital order as defined in kapdata
!
         JJ = 0
         DO K = 1, KAMAX
            IF (K > 1) JJ = NSHLII(K-1) + JJ
            DO J = 1, NSHLII(K)
               KK = NSHLPII(K,J)
               IF (NPY/=NP(KK) .OR. NAKY/=NAK(KK)) CYCLE
               JJJ = JJ + J
               PZII(JJJ) = PZY
               EII(JJJ) = EY
               NAKII(JJJ) = NAK(KK)
               NPII(JJJ) = NP(KK)
               NHII(JJJ) = NH(KK)
               GAMAII(JJJ) = GAMA(KK)
               CALL INTRPQI (PA, QA, MY, RA, JJJ, DNORM)
               IF (NTEST >= 100) WRITE (*, 301) NPII(JJJ), NHII(JJJ), EII(JJJ)&
                  , DNORM
               IF (NTEST > 1000) THEN
                  WRITE (*, *) 'PF              QF             RA'
                  DO KKK = 1, MFII(JJJ)
                     WRITE (*, *) PFII(KKK,JJJ), QFII(KKK,JJJ), RA(KKK)
                  END DO
               ENDIF
               NWIN = NWIN + 1
            END DO
         END DO
         CALL DALLOC (PA, 'PA', 'LODRWII')
         CALL DALLOC (QA, 'QA', 'LODRWII')
         CALL DALLOC (RA, 'RA', 'LODRWII')
         GO TO 3
      ENDIF
      IF (LDBPR(3)) WRITE (99, *) ' orbitals renormalised;'
!
!   Stop with an error message if all orbitals are not known
!
      IF (NWIN < NWII) THEN
         WRITE (6, *) 'LODRWF: All required orbitals not'
         WRITE (6, *) ' found.'
         IERR = 1
         GO TO 5
      ENDIF
!
      WRITE (6, *) ' ... load complete;'
!
    5 CONTINUE
      CLOSE(21)
      IF (NTEST >= 100) THEN
         WRITE (*, *) 'Sorted order should be the same as from kapdat'
         DO J = 1, NWII
            WRITE (*, 301) NPII(J), NHII(J), EII(J), DNORM
         END DO
         WRITE (*, *)
         WRITE (*, *) '*****************'
         WRITE (*, *) ' Leaving lodrwfi'
         WRITE (*, *) '*****************'
      ENDIF

      RETURN
!
  301 FORMAT(2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
      RETURN
!
      END SUBROUTINE LODRWFI
