!***********************************************************************
!
      SUBROUTINE GETOLDWT(NDEF, NCMIN, WT)
!
!   Interactively determines the weights for EOL calculation.
!   It's modified to always ask the question for the weight
!
!   Call(s) to: [LIB92]: GETYN.
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE iounit_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NDEF
      INTEGER , INTENT(IN) :: NCMIN
      REAL(DOUBLE)  :: WT(NCMIN)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, N159
      REAL(DOUBLE) :: SUMWGT
!-----------------------------------------------

! Standard weights: ncmin=1, OL calculation

      IF (NCMIN == 1) THEN
         WT(1) = -1.D0
         RETURN
      ENDIF

! Select a method to assign level weights for ncmin > 1 case

      WRITE (ISTDE, *) 'level weights (1 equal;  5 standard;  9 user)'

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      DO I = 1, 10
         READ (ISTDI, *) N159
         IF (N159==1 .OR. N159==5 .OR. N159==9) EXIT
         WRITE (ISTDE, *) 'Input not correct, do it again. tried=', I
      END DO

      IF (NDEF.EQ.0) THEN
        WRITE(734,*) n159,'! level weights'
      END IF

      IF (I > 10) STOP

!------------------------------------------------------------------

      SELECT CASE (N159)                         ! Equal weight
      CASE (1)
         WT = -2.D0
      CASE (5)                                   ! Standard weight
         WT = -1.D0
      CASE (9)                                   ! User-input weight

  123    CONTINUE
         WRITE (ISTDE, *) 'Enter the (relative) weights of the', NCMIN, &
            ' levels :'
         READ (ISTDI, *) (WT(I),I=1,NCMIN)

         SUMWGT = 0.D0
         DO I = 1, NCMIN
            IF (WT(I) <= 0.D0) THEN
               WRITE (ISTDE, *) 'Weights must exceed 0;'
               GO TO 123
            ELSE
               SUMWGT = SUMWGT + WT(I)
            ENDIF
         END DO
         SUMWGT = 1.D0/SUMWGT
         WT = SUMWGT*WT

      CASE DEFAULT
         WRITE (ISTDE, *) 'Impossible ! Because it was guarded'
         STOP
      END SELECT

      RETURN
      END SUBROUTINE GETOLDWT
