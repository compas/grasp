!***********************************************************************
!                                                                      *
      SUBROUTINE SET_CSF_ZFlist
!                                                                      *
!   Open, check, load data from and close the  .csl  file. This file   *
!   is always attached to stream 21.                                   *
!                                                                      *
!   Call(s) to:[RCI92]: OPENFL, PRSRSL, CONVRT, SET_CSF_NUMBER         *
!                                                                      *
!   Written by  G. Gaigalas                        Vilnius, May 2016   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
      USE IOUNIT_C
      USE ORB_C,               only: NW, NP, NAK
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I 
      USE prsrsl_I
      USE convrt_I
      USE Set_CSF_number_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: NCORE, IERR, IOS, I, J, LENTH, NAKJ, NPJ
      CHARACTER :: FILNAM*256, RECORD*15, FORM*11, STATUS*3
      CHARACTER :: RECORD_1*15, RECORD_2*15
      CHARACTER :: S_closed_1*181, S_closed_2*181
      CHARACTER :: S_orbitals_1*1070, S_orbitals_2*1070
!-----------------------------------------------
!
      print *,                                                          &
      "Give the full name of the list that contains the zero-order space"
      read(*,'(a)') FILNAM
      FORM = 'FORMATTED' 
      STATUS = 'OLD' 
      CALL OPENFL (21, FILNAM, FORM, STATUS, IERR) 
      IF (IERR == 1) THEN 
         WRITE (6, *) 'Error when opening', FILNAM 
         STOP  
      ENDIF 
      READ (21, '(1A15)', IOSTAT=IOS) RECORD_1
      IF (IOS/=0 .OR. RECORD_1(1:15)/='Core subshells:') THEN 
         WRITE (6, *) 'Not a Configuration Symmetry List File;' 
         CLOSE(21) 
      ENDIF 
!
      print *,                                                           &
      "Give the full name of the list that should be partitioned"
      read(*,'(a)') FILNAM
      WRITE (6, *) 'Loading Configuration Symmetry List File ...'
      FORM = 'FORMATTED' 
      STATUS = 'OLD' 
      CALL OPENFL (20, FILNAM, FORM, STATUS, IERR) 
      IF (IERR == 1) THEN 
         WRITE (6, *) 'Error when opening', FILNAM 
         STOP  
      ENDIF 
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
      IF (IOS/=0 .OR. RECORD_2(1:15)/='Core subshells:') THEN 
         WRITE (6, *) 'Not a Configuration Symmetry List File;' 
         CLOSE(20) 
      ENDIF 
!
      READ (21, '(A)') S_closed_1
      READ (20, '(A)') S_closed_2
      IF(S_closed_1 /= S_closed_2) then
         STOP "Diffeent close shells"
      end if
!
      FILNAM = 'rcsf.out' 
      FORM = 'FORMATTED' 
      STATUS = 'OLD' 
      CALL OPENFL (22, FILNAM, FORM, STATUS, IERR) 
      WRITE (22, '(1A15)') RECORD_1
      I = LEN_TRIM(S_closed_1)
      WRITE (22,'(A)') S_closed_1(1:I)
!
      READ (21, '(1A15)', IOSTAT=IOS) RECORD_1
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
      WRITE (22,'(1A15)') RECORD_1
!
      READ (21, '(A)') S_orbitals_1
      READ (20, '(A)') S_orbitals_2
!      I = LEN_TRIM(S_orbitals_1)
!      IF(S_orbitals_1(1:I) /= S_orbitals_2(1:I)) then
!         STOP "Diffeent order of Peel orbitals"
!      end if
      WRITE (22, '(A)') TRIM(S_orbitals_2)
!
      READ (21, '(1A7)', IOSTAT=IOS) RECORD_1
      READ (20, '(1A7)', IOSTAT=IOS) RECORD_2
      WRITE (22,'(1A7)') RECORD_1
!----------------------------------------------------------------------
      CALL Set_CSF_number
      REWIND(20)
      READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
!
!   Get the list of subshells
!
      NW = 0 
!
!   Read the list of core subshells; set up the arrays NP, NAK,
!   NKL, NKJ, NH for these subshells
!
      CALL PRSRSL (20, 1) 
      NCORE = NW 
!
!   Skip the peel subshell identification header; read the list of
!   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
!   these subshells
!
      READ (20, *) 
      CALL PRSRSL (20, 2) 
!
!   Ensure that the sets of core and peel subshell are disjoint
!
      DO J = NCORE + 1, NW 
         NPJ = NP(J) 
         NAKJ = NAK(J) 
         DO I = 1, NCORE 
            IF (NP(I)/=NPJ .OR. NAK(I)/=NAKJ) CYCLE  
            WRITE (ISTDE, *) 'SET_CSF_list: The lists of core and', &
               ' peel subshells must form disjoint sets.' 
            STOP  
         END DO 
      END DO 
!
!   Print the number of relativistic subshells
!
      IF (NW > 1) THEN 
         CALL CONVRT (NW, RECORD, LENTH) 
         WRITE (6, *) 'There are '//RECORD(1:LENTH)// &
                      ' relativistic subshells;' 
      ELSE 
         WRITE (6, *) 'There is 1 relativistic subshell;' 
      ENDIF 
      READ (20, *)
!
      RETURN  
      END SUBROUTINE SET_CSF_ZFlist
