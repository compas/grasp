!***********************************************************************
!
      SUBROUTINE GETALDWT(NCF, WT)
!
!   Interactively determines the weights.
!
!   Call(s) to: [LIB92]: ITJPO
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NCF
      REAL(DOUBLE)  :: WT(NCF)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, N159, NCMIN
      REAL(DOUBLE) :: SUMWGT, FTJPOI
!-----------------------------------------------

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

      IF (I > 10) STOP 'Must be running un-attended'

!------------------------------------------------------------------

      SELECT CASE (N159)                         ! Equal weight
      CASE (1)
         WT = 1.D0
         SUMWGT = DBLE(NCF)
      CASE (5)                                   ! Standard weight
         SUMWGT = 0.D0
         DO I = 1, NCF
            FTJPOI = DBLE(ITJPO(I))
            WT(I) = FTJPOI
            SUMWGT = SUMWGT + FTJPOI
         END DO
      CASE (9)                                   ! User-input weight

  123    CONTINUE
         WRITE (ISTDE, *) 'Enter the (relative) weights of the', NCF, &
            ' levels :'
         READ (ISTDI, *) (WT(I),I=1,NCMIN)

         SUMWGT = 0.D0
         DO I = 1, NCF
            IF (WT(I) <= 0.D0) THEN
               WRITE (ISTDE, *) 'Weights must exceed 0;'
               GO TO 123
            ELSE
               SUMWGT = SUMWGT + WT(I)
            ENDIF
         END DO

      CASE DEFAULT
         WRITE (ISTDE, *) 'Impossible ! Because it was guarded'
         STOP
      END SELECT

      SUMWGT = 1.D0/SUMWGT
      WT = SUMWGT*WT

      RETURN
      END SUBROUTINE GETALDWT
