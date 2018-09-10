!***********************************************************************
!                                                                      *
      SUBROUTINE LODRWF(IIERR) 
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!   Block version by Xinghong He          Last revision: 27 May 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:24:43   2/14/04  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNP
      USE memory_man
      USE DEBUG_C 
      USE DEF_C,           ONLY: C, Z 
      USE GRID_C 
      USE NPAR_C 
      USE ORB_C 
      USE WAVE_C,          ONLY: PZ, PF,QF 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE intrpq_I 
      USE orthsc_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: IIERR 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I, K, NWIN, IOS, NPY, NAKY, MY 
      integer :: ierr
      REAL(DOUBLE) :: CON, FKK, EY, PZY, DNORM 
      real(double), dimension(:), pointer :: PA, QA, RA
!-----------------------------------------------
!
!
!   Write entry message
!
      WRITE (6, *) 'Loading Radial WaveFunction File ...' 
!
!   Allocate storage to orbital arrays
!
      CALL ALLOC(PF, NNNP, NW, 'PF', 'LODRWF')
      CALL ALLOC(QF, NNNP, NW, 'QF', 'LODRWF')

!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      CON = Z/C 
      CON = CON*CON 
!
      DO J = 1, NW 
!
         PF(:N,J) = 0.0D00 
         QF(:N,J) = 0.0D00 
!
         E(J) = -1.0D00 
!
         K = ABS(NAK(J)) 
         IF (NPARM > 0) THEN 
            GAMA(J) = DBLE(K) 
         ELSE IF (NPARM == 0) THEN 
            FKK = DBLE(K*K) 
            IF (FKK >= CON) THEN 
               GAMA(J) = SQRT(FKK - CON) 
            ELSE 
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               STOP 'lodrwf: Inappropriate gamma' 
            ENDIF 
         ENDIF 
!
      END DO 
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      IF (LDBPR(3)) WRITE (99, 300) 
      NWIN = 0 
    3 CONTINUE 
      READ (23, IOSTAT=IOS) NPY, NAKY, EY, MY 
 
      IF (IOS == 0) THEN 
         CALL ALLOC (PA,MY, 'PA', 'LODRWF')
         CALL ALLOC (QA,MY, 'QA', 'LODRWF')
         CALL ALLOC (RA,MY, 'RA', 'LODRWF')
 
         READ (23) PZY, (PA(I),I=1,MY), (QA(I),I=1,MY) 
         READ (23) (RA(I),I=1,MY) 
 
         DO J = 1, NW 
            IF (.NOT.(E(J)<0.0D00 .AND. NPY==NP(J) .AND. NAKY==NAK(J))) CYCLE  
            PZ(J) = PZY 
            E(J) = EY 
            CALL INTRPQ(PA, QA, MY, RA, J, DNORM) 
            IF (LDBPR(3)) WRITE (99, 301) NP(J), NH(J), E(J), DNORM 
            NWIN = NWIN + 1 
         END DO 
         CALL DALLOC (PA, 'PA', 'LODRWF' )
         CALL DALLOC (QA, 'QA', 'LODRWF')
         CALL DALLOC (RA, 'RA', 'LODRWF')
         GO TO 3 
      ENDIF 
      IF (LDBPR(3)) WRITE (99, *) ' orbitals renormalised;' 
!
!   Stop with an error message if all orbitals are not known
!
      IF (NWIN < NW) THEN 
         IIERR = 1 
         GO TO 5 
      ENDIF 
!
!   Schmidt orthogonalise the orbitals
!
      CALL ORTHSC 
      IF (LDBPR(3)) WRITE (99, *) ' orbitals orthogonalised and renormalised;' 
!
      IIERR = 0 
    5 CONTINUE 
      RETURN  
!
  300 FORMAT(/,'From SUBROUTINE LODRWF:'/,' Orbital',8X,'Eigenvalue',19X,'Norm') 
  301 FORMAT(2X,I2,A2,4X,1P,1D22.15,4X,1D22.15) 
      RETURN  
!
      END SUBROUTINE LODRWF 
