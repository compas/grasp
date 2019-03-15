!***********************************************************************
!                                                                      *
      SUBROUTINE LODRWFI(NAME)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:25:11   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE memory_man
      USE def_C
      USE npar_C
      USE biorb_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER , INTENT(IN) :: NAME*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, K, NWIN, IOS, NPYII, NAKYII, MYII, JJ
      REAL(DOUBLE) :: CON, FKK, EYII, PZY
      REAL(DOUBLE), DIMENSION(:), pointer :: pa, qa, ra
      CHARACTER :: G92RWF*6
!-----------------------------------------------
!
!   Common relevant for the initial state
!
!   Write entry message
!
!      WRITE (6, *) 'Loading Radial WaveFunction File for initial state...'
!
!   Open the radial wave function file
!
      J = INDEX(NAME,' ')
      OPEN(UNIT=69, FILE=NAME(1:J-1)//'.bw', FORM='UNFORMATTED', STATUS='OLD', &
         POSITION='asis')
!
!   Allocate storage to orbital arrays
!
      CALL ALLOC (PFII, NNNP, NWII, 'PFII', 'LODRWFI')
      CALL ALLOC (QFII, NNNP, NWII, 'QFII', 'LODRWFI')
!
      CON = Z/C
      CON = CON*CON
!
      WRITE (*, *) 'NWII', NWII
      DO J = 1, NWII
         WRITE (*, *) NAKII(J), NPII(J), NHII(J)
         PFII(:NNNP,J) = 0.0D00
         QFII(:NNNP,J) = 0.0D00
!
         EII(J) = -1.0D00
!
         K = ABS(NAKII(J))
         IF (NPARM > 0) THEN
            GAMAII(J) = DBLE(K)
         ELSE IF (NPARM == 0) THEN
            FKK = DBLE(K*K)
            IF (FKK >= CON) THEN
               GAMAII(J) = SQRT(FKK - CON)
            ELSE
               WRITE (6, *) 'LODRWF: Imaginary gamma parameter'
               WRITE (6, *) ' for ', NPII(J), NHII(J), ' orbital; the'
               WRITE (6, *) ' point model for the nucleus'
               WRITE (6, *) ' is inappropriate for Z > ', C, '.'
               STOP
            ENDIF
         ENDIF
!
      END DO
!
!   Read orbital information from Read Orbitals File;
!
      NWIN = 0
      READ (69, IOSTAT=IOS) G92RWF
      IF (IOS/=0 .OR. G92RWF/='G92RWF') THEN
         WRITE (6, *) 'This is not a Radial WaveFunction File;'
         CLOSE(69)
      ENDIF

    3 CONTINUE
      READ (69, IOSTAT=IOS) NPYII, NAKYII, EYII, MYII
      IF (IOS == 0) THEN
         CALL ALLOC (PA, MYII, 'PA', 'LODRWFI')
         CALL ALLOC (QA, MYII, 'QA', 'LODRWFI')
         CALL ALLOC (RA, MYII, 'RA', 'LODRWFI')
         READ (69) PZY, (PA(I),I=1,MYII), (QA(I),I=1,MYII)
         READ (69) (RA(I),I=1,MYII)

         DO J = 1, NWII
            IF (.NOT.(EII(J)<0.0D00 .AND. NPYII==NPII(J) .AND. NAKYII==NAKII(J)&
               )) CYCLE
            PZII(J) = PZY
            EII(J) = EYII
            MFII(J) = MYII
            DO JJ = 1, MFII(J)
               PFII(JJ,J) = PA(JJ)
               QFII(JJ,J) = QA(JJ)
            END DO
            NWIN = NWIN + 1
         END DO
         CALL DALLOC (PA, 'PA', 'LODRWFI')
         CALL DALLOC (QA, 'QA', 'LODRWFI')
         CALL DALLOC (RA, 'RA', 'LODRWFI')
         GO TO 3
      ENDIF
!
!   Stop with an error message if all orbitals are not known
!
      IF (NWIN < NWII) THEN
         WRITE (6, *) 'LODRWF: All required orbitals not'
         WRITE (6, *) ' found.'
         STOP
      ENDIF
!
!      WRITE (6, *) ' ... load complete;'
!
      CLOSE(69)

      RETURN
      END SUBROUTINE LODRWFI
