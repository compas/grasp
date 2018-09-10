!***********************************************************************
!                                                                      *
      SUBROUTINE KAPDATA(NTESTG, NCORE1, NCORE2) 
!                                                                      *
!   This subroutine determines the number of kappa quantum numbers     *
!   KAMAX together with the number of orbitals of each kappa.          *
!                                                                      *
!   This subroutine also checks if the orbital ordering is normal      *
!   or reversed. If reversed then JA and JB have to be permuted        *
!   In the subroutine ti1tv                                            *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:26:44   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE biorb_C,         ONLY: nwii,nakii,npii,nwff,nakff,npff
      USE sbdat_C,         ONLY: NLMAX,NINII,NINFF,NSHLII,NSHLFF,      &
                                 IKAPPA,NAKINVII,NAKINVFF,NLMAX,       &
                                 NSHLPII,NSHLPFF,NSHLPPII,NSHLPPFF,KAMAX
      USE orbord_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NTESTG 
      INTEGER, INTENT(IN) :: NCORE1 
      INTEGER, INTENT(IN) :: NCORE2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1 
      INTEGER, DIMENSION(2*NNNW) :: ISORT 
      INTEGER :: NTESTL, NTEST, I, J, K, ITAL, NREF 
      LOGICAL :: KLAR 
!-----------------------------------------------
!
!
!
      NTESTL = 0 
      NTEST = MAX0(NTESTL,NTESTG) 
 
      ISORT(:2*NNNW) = 100 
 
      IKAPPA(:NLMAX) = 0 
      NSHLII(:NLMAX) = 0 
      NSHLFF(:NLMAX) = 0 
      NINII(:NLMAX) = 0 
      NINFF(:NLMAX) = 0 
 
      NAKINVII(:NNNW) = 0 
      NAKINVFF(:NNNW) = 0 
 
      NSHLPII(:NLMAX,:NLMAX) = 0 
      NSHLPFF(:NLMAX,:NLMAX) = 0 
 
      NSHLPPII(:NLMAX,:NNNW) = 0 
      NSHLPPFF(:NLMAX,:NNNW) = 0 
!
!   Sort the kappa quantum numbers
!
      DO K = 1, NWII + NWFF 
         IF (K <= NWII) THEN 
            ITAL = NAKII(K) 
         ELSE 
            ITAL = NAKFF(K-NWII) 
         ENDIF 
         I = K - 1 
         KLAR = .FALSE. 
   12    CONTINUE 
         IF (I>0 .AND. .NOT.KLAR) THEN 
            IF (ITAL <= ISORT(I)) THEN 
               ISORT(I+1) = ISORT(I) 
               I = I - 1 
            ELSE 
               KLAR = .TRUE. 
            ENDIF 
            GO TO 12 
         ENDIF 
         ISORT(I+1) = ITAL 
      END DO 
!
!   Determine the unique set of kappa IKAPPA
!
      KAMAX = 1 
      IKAPPA(1) = ISORT(1) 
      DO K = 1, 2*NNNW - 1 
         IF (ISORT(K)==ISORT(K+1) .OR. ISORT(K+1)>=100) CYCLE  
         KAMAX = KAMAX + 1 
         IKAPPA(KAMAX) = ISORT(K+1) 
      END DO 
!
!  Make a connection between each kappa and a number in the
!  range [1,KAMAX] as to know on which file to dump the data
!  Determine the number of shells NSHLII for each I in the
!  range [1,KAMAX]
!
      IF (NTEST >= 10) THEN 
         WRITE (*, *) '******************' 
         WRITE (*, *) ' Entering kapdata' 
         WRITE (*, *) '******************' 
         WRITE (*, *) 
         WRITE (*, *) 'There are', NWII, 'orbitals in the initial state' 
         WRITE (*, *) 'with the following n and kappa quantum numbers' 
      ENDIF 
 
      DO J = 1, NWII 
         IF (NTEST >= 10) WRITE (*, *) 'orbital number', J, 'n and kappa', NPII&
            (J), NAKII(J) 
         IF (J <= NCORE1) THEN 
            DO I = 1, KAMAX 
               IF (IKAPPA(I) /= NAKII(J)) CYCLE  
               NAKINVII(J) = I 
               NINII(I) = NINII(I) + 1 
               NSHLII(I) = NSHLII(I) + 1 
               NSHLPII(I,NSHLII(I)) = J 
               NSHLPPII(I,J) = NSHLII(I) 
            END DO 
         ELSE 
            DO I = 1, KAMAX 
               IF (IKAPPA(I) /= NAKII(J)) CYCLE  
               NAKINVII(J) = I 
               NSHLII(I) = NSHLII(I) + 1 
               NSHLPII(I,NSHLII(I)) = J 
               NSHLPPII(I,J) = NSHLII(I) 
            END DO 
         ENDIF 
      END DO 
 
      IF (NTEST >= 10) THEN 
         WRITE (*, *) 'There are', NWFF, 'orbitals in the final state' 
         WRITE (*, *) 'with the following n and kappa quantum numbers' 
      ENDIF 
      DO J = 1, NWFF 
         IF (NTEST >= 10) WRITE (*, *) 'orbital number', J, 'n and kappa=', &
            NPFF(J), NAKFF(J) 
         IF (J <= NCORE2) THEN 
            DO I = 1, KAMAX 
               IF (IKAPPA(I) /= NAKFF(J)) CYCLE  
               NAKINVFF(J) = I 
               NINFF(I) = NINFF(I) + 1 
               NSHLFF(I) = NSHLFF(I) + 1 
               NSHLPFF(I,NSHLFF(I)) = J 
               NSHLPPFF(I,J) = NSHLFF(I) 
            END DO 
         ELSE 
            DO I = 1, KAMAX 
               IF (IKAPPA(I) /= NAKFF(J)) CYCLE  
               NAKINVFF(J) = I 
               NSHLFF(I) = NSHLFF(I) + 1 
               NSHLPFF(I,NSHLFF(I)) = J 
               NSHLPPFF(I,J) = NSHLFF(I) 
            END DO 
         ENDIF 
      END DO 
 
      IF (NTEST >= 10) THEN 
         WRITE (*, *) 'Total number of different kappa', KAMAX 
         DO I = 1, KAMAX 
            WRITE (*, *) 'L=', I, 'corresponds to kappa=', IKAPPA(I) 
            WRITE (*, *) 'nr of init.  orb. with this kappa=', NSHLII(I) 
            WRITE (*, *) 'nr of final. orb. with this kappa=', NSHLFF(I) 
         END DO 
         DO I = 1, KAMAX 
            WRITE (*, *) 'Position in initial state list' 
            DO J = 1, NSHLII(I) 
               WRITE (*, *) 'L=', I, 'orb. nr', J, ',position', NSHLPII(I,J) 
            END DO 
            WRITE (*, *) 'Position in final   state list' 
            DO J = 1, NSHLFF(I) 
               WRITE (*, *) 'L=', I, 'orb. nr', J, ',position', NSHLPFF(I,J) 
            END DO 
            WRITE (*, *) 'Relative positions for initial state orbitals' 
            DO J = 1, NWII 
               IF (NSHLPPII(I,J) == 0) CYCLE  
               WRITE (*, *) 'Orbital', J, 'is nr', NSHLPPII(I,J), 'with kappa'&
                  , IKAPPA(I) 
            END DO 
            WRITE (*, *) 'Relative positions for final state orbitals' 
            DO J = 1, NWFF 
               IF (NSHLPPFF(I,J) == 0) CYCLE  
               WRITE (*, *) 'Orbital', J, 'is nr', NSHLPPFF(I,J), 'with kappa'&
                  , IKAPPA(I) 
            END DO 
         END DO 
 
      ENDIF 
!
!  Check if the orbital ordering is normal or reversed.
!
!  For normal ordering
!
!    NPII(NSHLPII(I,1)) < NPII(NSHLPII(I,2)) < ... < NPII(NSHLPII(I,NSHLII(I))
!
!  For reversed ordering
!
!    NPII(NSHLPII(I,1)) > NPII(NSHLPII(I,2)) > ... > NPII(NSHLPII(I,NSHLII(I))
!
      NORDII = 0 
      DO I = 1, KAMAX 
         NREF = 0 
         DO J = 1 + NINII(I), NSHLII(I) 
            IF (NPII(NSHLPII(I,J)) < NREF) NORDII = 1 
            NREF = NPII(NSHLPII(I,J)) 
         END DO 
      END DO 
 
      NORDFF = 0 
      DO I = 1, KAMAX 
         NREF = 0 
         DO J = 1 + NINFF(I), NSHLFF(I) 
            IF (NPFF(NSHLPFF(I,J)) < NREF) NORDFF = 1 
            NREF = NPFF(NSHLPFF(I,J)) 
         END DO 
      END DO 
 
      IF (NORDII /= NORDFF) THEN 
         WRITE (*, *) ' Orbital order of the initial and final states' 
         WRITE (*, *) ' should be the same. STOP' 
         STOP  
      ENDIF 
!
!   If not normal order check if reversed order
!
      IF (NORDII == 1) THEN 
         DO I = 1, KAMAX 
            NREF = 0 
            DO J = NSHLII(I), 1 + NINII(I), -1 
               IF (NPII(NSHLPII(I,J)) < NREF) NORDII = 2 
               NREF = NPII(NSHLPII(I,J)) 
            END DO 
         END DO 
 
         DO I = 1, KAMAX 
            NREF = 0 
            DO J = NSHLFF(I), 1 + NINFF(I), -1 
               IF (NPFF(NSHLPFF(I,J)) < NREF) NORDFF = 2 
               NREF = NPFF(NSHLPFF(I,J)) 
            END DO 
         END DO 
      ENDIF 
 
      IF (NORDII==2 .OR. NORDFF==2) THEN 
         WRITE (*, *) ' The orbital order is neither normal or reversed' 
         WRITE (*, *) ' STOP' 
         STOP  
      ENDIF 
 
!w      write(*,*) 'Give nordii'
!w      read(*,*) nordii
!w      nordff = nordii
 
      IF (NORDII == 0) THEN 
         WRITE (*, *) ' Normal orbital ordering' 
      ELSE 
         WRITE (*, *) ' Reverse orbital ordering' 
      ENDIF 
 
      WRITE (*, *) '*****************' 
      WRITE (*, *) ' Leaving kapdata' 
      WRITE (*, *) '*****************' 
      WRITE (*, *) 
 
      RETURN  
      END SUBROUTINE KAPDATA 
