!***********************************************************************
!                                                                      *
      SUBROUTINE NCHARG
!                                                                      *
!   This routine evaluates the nuclear charge density, and stores it   *
!   in the  common  array  ZDIST .                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ES.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE def_C,           ONLY: pi,  z, precis
      USE grid_C
      USE npar_C
      USE ncdist_C
      USE tatb_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE es_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I
      REAL(DOUBLE) :: C, A, CBA, PI2, ABC, ABC2, ABC3, S2MCBA, S3MCBA, EN,&
         ZNORM, EXTRM, ZDISTI
      LOGICAL :: FORM1, FORM2
!-----------------------------------------------
!
!
!   Initialize array to zero
!
      ZDIST(:N) = 0.0D00
!
!   Fermi charge distribution
!
      IF (NPARM == 2) THEN
         C = PARM(1)
         A = PARM(2)
         CBA = C/A
         PI2 = PI*PI
         ABC = A/C
         ABC2 = ABC*ABC
         ABC3 = ABC2*ABC
         CALL ES ((-CBA), S2MCBA, S3MCBA)
         EN = 1.0D00 + PI2*ABC2 - 6.0D00*ABC3*S3MCBA
         ZNORM = 3.0D00*Z/(4.0D00*PI*EN*C**3)
         FORM1 = .TRUE.
         FORM2 = .FALSE.
         DO I = 1, N
            IF (FORM1) THEN
               EXTRM = DEXP((R(I)-C)/A)
               ZDIST(I) = ZNORM/(1.0D00 + EXTRM)
               IF (1.0D00/EXTRM <= PRECIS) THEN
                  FORM1 = .FALSE.
                  FORM2 = .TRUE.
               ENDIF
            ELSE IF (FORM2) THEN
               ZDISTI = ZNORM*DEXP((-(R(I)-C)/A))
               IF (DABS(ZDISTI) > 0.0D00) THEN
                  ZDIST(I) = ZDISTI
               ELSE
                  MTP = I
                  EXIT
               ENDIF
            ENDIF
         END DO
      ENDIF
!
      RETURN
!
      END SUBROUTINE NCHARG
