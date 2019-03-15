!***********************************************************************
!                                                                      *
      SUBROUTINE FNAME(NAME)
!                                                                      *
!   Determines the name of the initial and final states                *
!   In addition this subroutine determines which J symmetries          *
!   that are to be transformed                                         *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:19:37   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   C O M M O N   B l o c k s
!-----------------------------------------------
      USE jqjc_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I
      LOGICAL :: YES
!-----------------------------------------------
!
!
!   Obtain the names of the initial and final state files
!
    1 CONTINUE
      WRITE (6, *) ' Name of the Initial state'
      READ (*, '(A)') NAME(1)

      WRITE (6, *) ' Name of the Final state'
      READ (*, '(A)') NAME(2)
!
      J = INDEX(NAME(1),' ')
      IF (J == 1) THEN
         WRITE (6, *) ' Names may not start with blanks'
         GO TO 1
      ENDIF
!
      J = INDEX(NAME(2),' ')
      IF (J == 1) THEN
         WRITE (6, *) ' Names may not start with blanks'
         GO TO 1
      ENDIF
! Per april 2007
! Check if the initial and final states are identical.

      IF (TRIM(NAME(1)) == TRIM(NAME(2))) THEN
         PRINT *
         PRINT *, ' Initial and final states are identical and there is'
         PRINT *, ' no need for the biorthogonal transformation. Just  '
         PRINT *, ' copy name.w to name.bw and name.(c)m to name.(c)bm '
         PRINT *, ' and run rtransition.                               '
         PRINT *


         PRINT *, ' Do you want to continue anyway ? '
         YES = GETYN ()
         IF (YES.EQV..FALSE.) STOP
      END IF

! end Per 2007
      WRITE (6, *) ' Transformation of all J symmetries?'
      YES = GETYN()
      IF (YES) THEN
         NTRANS = 0
      ELSE
         NTRANS = 1
         WRITE (6, *) ' Number of initial state J symmetries to be transformed'
         READ (*, *) JQJ1
         WRITE (6, *) ' Give the J symmetries in the form 2*J'
         READ (*, *) (ITJQJ1(I),I=1,JQJ1)
         ITJQJ1(:JQJ1) = ITJQJ1(:JQJ1) + 1
         WRITE (6, *) ' Number of final state J symmetries to be transformed'
         READ (*, *) JQJ2
         WRITE (6, *) ' Give the J symmetries in the form 2*J'
         READ (*, *) (ITJQJ2(I),I=1,JQJ2)
         ITJQJ2(:JQJ2) = ITJQJ2(:JQJ2) + 1
      ENDIF

      RETURN
      END SUBROUTINE FNAME
