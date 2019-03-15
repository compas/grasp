!***********************************************************************
!                                                                      *
      SUBROUTINE MERG12(NAME, NCORER, NCORE)
!                                                                      *
!   This subroutines merges the initial and final state lists          *
!   Observ that there may doublets in this list if the initial         *
!   and final states have the same parity                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE def_C
      USE orb_C
      USE blk_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE iq_I
      USE iqr_I
      USE ispar_I
      USE isparr_I
      USE itjpo_I
      USE itjpor_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NCORER
      INTEGER , INTENT(IN) :: NCORE
      CHARACTER , INTENT(IN) :: NAME(2)*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1
      INTEGER, DIMENSION(0:NNNW + 1) :: NEWNP
      INTEGER, DIMENSION(NNNW) :: NEWNAK, ICON, ICON1
      INTEGER :: I, J, ICOMP, NEW, IPI1, IPI2, K, NEW1, NEW2, M, NLINE
      CHARACTER :: LINE*500
      CHARACTER, DIMENSION(NNNW) :: NEWNH*2
      CHARACTER :: NEW3*2
!-----------------------------------------------
!
!
! NCFI(I): the end position of the Ith block for the initial states in the globle CSF list
! NCFF(I): the end position of the Ith block for the final states in the globle CSF list

!
      OPEN(UNIT=21, FILE='SLASK', FORM='FORMATTED', STATUS='UNKNOWN', POSITION=&
         'asis')
!
!   The same number of electrons must appear in both lists
!
      IF (NELECR /= NELEC) THEN
         WRITE (6, *) 'The number of electrons is not equal in the'
         WRITE (6, *) ' first and second GRASP92 Configuration'
         WRITE (6, *) ' symmetry list files.'
         STOP
      ENDIF
!
!   THe core orbitals must be the same
!
      IF (NCORE /= NCORER) THEN
         WRITE (6, *) 'The number of core orbitals must be the same'
         STOP
      ENDIF

      DO I = 1, NCORE
         IF (NP(I)==NPR(I) .AND. NAK(I)==NAKR(I)) CYCLE
         WRITE (6, *) 'The core orbitals must be the same'
         STOP
      END DO

      ICON(:NW) = 0
!
!   For each orbital in the initial list check if there is a corresponding
!   orbital in the final list. If so give the number of the corresponding
!   orbital
!
      DO I = 1, NW
         DO J = 1, NWR
            IF (NP(I)/=NPR(J) .OR. NAK(I)/=NAKR(J)) CYCLE
            ICON(I) = J
         END DO
      END DO
!
!   Check if the ordering of the initial and final state orbitals
!   is consistent. The condition for this is that ICON is increasing
!
      J = 0
      DO I = 1, NW
         IF (ICON(I) == 0) CYCLE
         J = J + 1
         ICON1(J) = ICON(I)
      END DO

      ICOMP = ICON1(1)
      DO I = 2, J
         IF (ICON1(I) < ICOMP) THEN
            WRITE (*, *) ' In merg12: ordering of the initial and final'
            WRITE (*, *) ' state orbitals is inconsistent. STOP'
            STOP
         ELSE
            ICOMP = ICON1(I)
         ENDIF
      END DO
!
!   Determine a common orbital set for the initial and final state.
!   The common set must be such that the order of both the initial and
!   final state sets are preserved.
!
      DO I = 1, NWR
         NEWNP(I) = NPR(I)
         NEWNAK(I) = NAKR(I)
         NEWNH(I) = NHR(I)
         WRITE (*, *) NEWNP(I), NEWNH(I)
      END DO
!
!   Add the initial state orbitals at the end
!
      NEW = NWR
      DO I = 1, NW
         IF (ICON(I) /= 0) CYCLE
         NEW = NEW + 1
         NEWNP(NEW) = NP(I)
         NEWNAK(NEW) = NAK(I)
         NEWNH(NEW) = NH(I)
         WRITE (*, *) NEWNP(NEW), NEWNH(NEW), NEW
      END DO
!
!   Now sort in the orbitals at the end in the right position
!
      L193: DO I = NWR + 1, NEW
!
!   Position in initial state list
!
         DO J = 1, NW
            IF (NEWNP(I)/=NP(J) .OR. NEWNAK(I)/=NAK(J)) CYCLE
            IPI1 = J
         END DO
         WRITE (*, *) 'i,ipi1', I, IPI1

         DO J = 1, NWR
            IPI2 = 0
            DO K = 1, NW
               IF (NEWNP(J)/=NP(K) .OR. NEWNAK(J)/=NAK(K)) CYCLE
               IPI2 = K
            END DO
            WRITE (*, *) 'j,ipi2', J, IPI2
            IF (IPI2 == 0) CYCLE
            IF (IPI1 >= IPI2) CYCLE
            NEW1 = NEWNP(I)
            NEW2 = NEWNAK(I)
            NEW3 = NEWNH(I)
            NEWNP(I:1+J:(-1)) = NEWNP(I-1:J:(-1))
            NEWNAK(I:1+J:(-1)) = NEWNAK(I-1:J:(-1))
            NEWNH(I:1+J:(-1)) = NEWNH(I-1:J:(-1))
            NEWNP(J) = NEW1
            NEWNAK(J) = NEW2
            NEWNH(J) = NEW3
            DO M = 1, NEW
               WRITE (*, *) NEWNP(M), NEWNH(M)
            END DO
            CYCLE  L193
         END DO
      END DO L193

      NW = NEW
      NP(:NW) = NEWNP(1:NW)
      NAK(:NW) = NEWNAK(:NW)
      NH(:NW) = NEWNH(:NW)
!
!   Determine NKL and NKJ
!
      DO I = 1, NW
         NKJ(I) = 2*ABS(NAK(I)) - 1
         IF (NAK(I) > 0) THEN
            NKL(I) = (NKJ(I)+1)/2
         ELSE
            NKL(I) = (NKJ(I)-1)/2
         ENDIF
      END DO
!
!   Determine the common core subshells; write out the list;
!   determine the pell subshells; write out the list; these
!   data form the first part of the header of the .csl file;
!   one additional line forms the remainder of the header of
!   the  .csl file
!
      WRITE (21, '(A)') 'Core subshells:'
      WRITE (21, 301) (NP(I),NH(I),I=1,NCORE)
      WRITE (21, '(A)') 'Peel subshells:'
      WRITE (21, 301) (NP(I),NH(I),I=NCORE + 1,NW)
      WRITE (21, '(A)') 'CSF(s):'
!
!   Now write out all CSFs in the initial and final state list
!
      J = INDEX(NAME(1),' ')
      OPEN(UNIT=23, FILE=NAME(1)(1:J-1)//'.c', FORM='FORMATTED', STATUS='OLD', &
         POSITION='asis')

      DO I = 1, 5
         READ (23, '(A)') LINE
      END DO

      NBLOCKI = 0
      NLINE = 0
    5 CONTINUE
      READ (23, '(A)', END=98) LINE
      IF (LINE(1:2) == ' *') THEN
         NBLOCKI = NBLOCKI + 1
         NCFI(NBLOCKI) = NLINE/3
      ELSE
         NLINE = NLINE + 1
      ENDIF
      K = 500
   10 CONTINUE
      IF (LINE(K:K) == ' ') THEN
         K = K - 1
         IF (K > 1) GO TO 10
      ENDIF
      WRITE (21, '(A)') LINE(1:K)
      GO TO 5
   98 CONTINUE
      NBLOCKI = NBLOCKI + 1
      NCFI(NBLOCKI) = NLINE/3
      CLOSE(23)

! zou
!     if(NAME(2).EQ.NAME(1)) return
! zou
      J = INDEX(NAME(2),' ')
      OPEN(UNIT=23, FILE=NAME(2)(1:J-1)//'.c', FORM='FORMATTED', STATUS='OLD', &
         POSITION='asis')

      DO I = 1, 5
         READ (23, '(A)') LINE
      END DO

      NBLOCKF = 0
      NLINE = 0
   15 CONTINUE
      READ (23, '(A)', END=99) LINE
      IF (LINE(1:2) == ' *') THEN
         NBLOCKF = NBLOCKF + 1
         NCFF(NBLOCKF) = NLINE/3
      ELSE
         NLINE = NLINE + 1
      ENDIF
      K = 500
   20 CONTINUE
      IF (LINE(K:K) == ' ') THEN
         K = K - 1
         IF (K > 1) GO TO 20
      ENDIF
      WRITE (21, '(A)') LINE(1:K)
      GO TO 15
   99 CONTINUE
      NBLOCKF = NBLOCKF + 1
      NCFF(NBLOCKF) = NLINE/3
      CLOSE(23)

      CLOSE(21)
!

      WRITE (6, *) NBLOCKI
      WRITE (6, *) (NCFI(I),I=1,NBLOCKI)
      WRITE (6, *) NBLOCKF
      WRITE (6, *) (NCFF(I),I=1,NBLOCKF)
  301 FORMAT(120(1X,1I2,1A2))
!
      RETURN
      END SUBROUTINE MERG12
