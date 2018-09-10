!***********************************************************************
!                                                                      *
      SUBROUTINE LODRWFF(NAME) 
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
      CHARACTER , INTENT(IN) :: NAME*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, K, NWIN, IOS, NPYFF, NAKYFF, MYFF, JJ 
      REAL(DOUBLE) :: CON, FKK, EYFF, PZY 
      REAL(DOUBLE), DIMENSION(:), pointer :: pa, qa, ra
      CHARACTER :: G92RWF*6 
!-----------------------------------------------
!
!   Common relevant for the final state
!
!
!   Write entry message
!
      WRITE (6, *) 'Loading Radial WaveFunction File for final state...' 
!
!   Open the radial wave function file
!
      J = INDEX(NAME,' ') 
      OPEN(UNIT=69, FILE=NAME(1:J-1)//'.bw', FORM='UNFORMATTED', STATUS='OLD', &
         POSITION='asis') 
!
!   Allocate storage to orbital arrays
!
      CALL ALLOC (PFFF, NNNP, NWFF, 'PFFF', 'LODRWFF') 
      CALL ALLOC (QFFF, NNNP, NWFF, 'QFFF', 'LODRWFF') 
!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      CON = Z/C 
      CON = CON*CON 
!
      DO J = 1, NWFF 
         PFFF(:NNNP,J) = 0.0D00 
         QFFF(:NNNP,J) = 0.0D00 
!
         EFF(J) = -1.0D00 
!
         K = ABS(NAKFF(J)) 
         IF (NPARM > 0) THEN 
            GAMAFF(J) = DBLE(K) 
         ELSE IF (NPARM == 0) THEN 
            FKK = DBLE(K*K) 
            IF (FKK >= CON) THEN 
               GAMAFF(J) = SQRT(FKK - CON) 
            ELSE 
               WRITE (6, *) 'LODRWF: Imaginary gamma parameter' 
               WRITE (6, *) ' for ', NPFF(J), NHFF(J), ' orbital; the' 
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
      READ (69, IOSTAT=IOS) NPYFF, NAKYFF, EYFF, MYFF 
      IF (IOS == 0) THEN 
         CALL ALLOC (PA, MYFF, 'PA', 'LODRWFF' ) 
         CALL ALLOC (QA, MYFF, 'QA', 'LODRWFF' ) 
         CALL ALLOC (RA, MYFF, 'RA', 'LODRWFF' ) 
         READ (69) PZY, (PA(I),I=1,MYFF), (QA(I),I=1,MYFF) 
         READ (69) (RA(I),I=1,MYFF) 
 
         DO J = 1, NWFF 
            IF (.NOT.(EFF(J)<0.0D00 .AND. NPYFF==NPFF(J) .AND. NAKYFF==NAKFF(J)&
               )) CYCLE  
            PZFF(J) = PZY 
            EFF(J) = EYFF 
            MFFF(J) = MYFF 
            DO JJ = 1, MFFF(J) 
               PFFF(JJ,J) = PA(JJ) 
               QFFF(JJ,J) = QA(JJ) 
            END DO 
            NWIN = NWIN + 1 
         END DO 
         CALL DALLOC (PA, 'PA', 'LODRWFF') 
         CALL DALLOC (QA, 'QA', 'LODRWFF') 
         CALL DALLOC (RA, 'RA', 'LODRWFF') 
         GO TO 3 
      ENDIF 
!
!   Stop with an error message if all orbitals are not known
!
      IF (NWIN < NWFF) THEN 
         WRITE (6, *) 'LODRWF: All required orbitals not' 
         WRITE (6, *) ' found.' 
         STOP  
      ENDIF 
!
      WRITE (6, *) ' ... load complete;' 
!
      CLOSE(69) 
 
      RETURN  
      END SUBROUTINE LODRWFF 
