

!***********************************************************************
!                                                                      *
      SUBROUTINE GENRWF
!                                                                      *
!   Controls the computation of the subshell radial wavefunctions.     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETRSL, GETYN.                         *
!               [ERWF]: FRMHYD, FRMRWF, FRMTFP, PRTREM, SUMMRY,        *
!                       TFPOT.                                         *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Updated by Xinghong He                Last revision: 11 Nov 1997   *
!    for menu-driven and less questions
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE COUN_C
      USE DEF_C, ONLY: Z, C, NNNP
      USE DEFAULT_C
      USE GRID_C
      USE LEFT_C, ONLY: SET
      USE NPAR_C
      USE ORB_C
      USE WAVE_C, ONLY: PF, QF
      USE WHFROM_C, ONLY: SOURCE
      USE IOUNIT_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE tfpot_I
      USE prtrem_I
      USE getrsl_I
      USE frmrwf_I
      USE frmtfp_I
      USE frmhyd_I
      USE summry_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(NNNW) :: INDEX
      INTEGER :: J, I, K, NRADIAL, NSUBS, LOC
      REAL(DOUBLE) :: CON, FKK
      LOGICAL :: ALL, YES, MODIFY=.FALSE., EXISTED
      CHARACTER :: INFILE*128
!-----------------------------------------------

!
!   Set the threshold for node counting
!
      THRESH = 0.05
!
!   Set up the Thomas-Fermi potential
!
      CALL TFPOT
!
!   Allocate storage to the arrays that store the subshell
!   radial wavefunction arrays; initialise these and all
!   associated arrays
!
      CALL ALLOC (PF,NNNP,NW,'PF', 'GENRWF')
      CALL ALLOC (QF,NNNP,NW,'QF', 'GENRWF')

      CON = Z/C
      CON = CON*CON

      DO J = 1, NW
         SET(J) = .FALSE.
         SOURCE(J) = ' '
         PF(:N,J) = 0.D0
         QF(:N,J) = 0.D0

         K = ABS(NAK(J))
         IF (NPARM > 0) THEN
            GAMA(J) = DBLE(K)
         ELSE IF (NPARM == 0) THEN
            FKK = DBLE(K*K)
            IF (FKK >= CON) THEN
               GAMA(J) = SQRT(FKK - CON)
            ELSE
               WRITE (ISTDE, *) 'LODRWF: Imaginary gamma parameter ', 'for ', &
                  NP(J), NH(J), ' orbital;'
               WRITE (ISTDE, *) 'the point model for the nucleus ', &
                  'is inappropriate for Z > ', C, '.'
               STOP
            ENDIF
         ENDIF
      END DO
!
!   Write out the complete list of subshell radial wave functions
!
      CALL PRTREM (ALL)
!
! Direct to read radial functions till finish
!
  123 CONTINUE
      IF (.NOT.ALL) THEN
  234    CONTINUE
         WRITE (ISTDE, *)
         WRITE (ISTDE, *) 'Read subshell radial wavefunctions. ', &
            'Choose one below'
         WRITE (ISTDE, *) '    1 -- GRASP92 File'
         WRITE (ISTDE, *) '    2 -- Thomas-Fermi'
         WRITE (ISTDE, *) '    3 -- Screened Hydrogenic'
         WRITE (ISTDE, *) '    4 -- Screened Hydrogenic [custom Z]'

         READ (ISTDI, *) NRADIAL
         IF (NRADIAL<1 .OR. NRADIAL>4) THEN
            WRITE (ISTDE, *) NRADIAL, 'is not a valid choice, redo'
            GO TO 234
         ENDIF

         IF (NRADIAL == 1) THEN
  345       CONTINUE
            WRITE (ISTDE, *) 'Enter the file name (Null then "rwfn.out")'
            READ (ISTDI, '(A)') INFILE
            IF (LEN_TRIM(INFILE) == 0) INFILE = 'rwfn.out'

            INQUIRE(FILE=INFILE, EXIST=EXISTED)
            IF (.NOT.EXISTED) THEN
               WRITE (ISTDE, *) ' File "', INFILE(1:LEN_TRIM(INFILE)), &
                  '" does not exist, redo'
               GO TO 345
            ENDIF
         ENDIF

         WRITE (ISTDE, *) 'Enter the list of relativistic subshells:'
         OPEN(UNIT=734, FILE='tmp_734', STATUS='NEW')
         CALL GETRSL (INDEX, NSUBS)
         CLOSE(734, STATUS='DELETE')

         IF (NRADIAL == 1) THEN
            CALL FRMRWF (INDEX, NSUBS, INFILE)
         ELSE IF (NRADIAL == 2) THEN
            CALL FRMTFP (INDEX, NSUBS)
         ELSE IF (NRADIAL == 3) THEN
            CALL FRMHYD (INDEX, NSUBS, MODIFY,NRADIAL)
         ELSE IF (NRADIAL == 4) THEN
            CALL FRMHYD (INDEX, NSUBS, MODIFY,NRADIAL)
         ENDIF

         CALL PRTREM (ALL)
         IF (.NOT.ALL) GO TO 234

      ENDIF
!
! All read. Let know, and allow modifying if non default
!
      WRITE (ISTDE, *) 'All required subshell radial wavefunctions ', &
         ' have been estimated:'
      CALL SUMMRY (ISTDE)

      IF (NDEF == 0) THEN
         MODIFY = .FALSE.
      ELSE
         WRITE (ISTDE, *) 'Revise any of these estimates?'
         MODIFY = GETYN()
      ENDIF

      IF (MODIFY) THEN
  456    CONTINUE
         WRITE (ISTDE, *) 'Enter the list of subshells whose radial ', &
            'wavefunctions are to be revised:'
         OPEN(UNIT=734, FILE='tmp_734', STATUS='NEW')
         CALL GETRSL (INDEX, NSUBS)
         CLOSE(734, STATUS='DELETE')
         IF (NSUBS == 0) GO TO 456
         DO J = 1, NSUBS
            LOC = INDEX(J)
            SET(LOC) = .FALSE.
            PF(:N,LOC) = 0.D0
            QF(:N,LOC) = 0.D0
            SOURCE(LOC) = ' '
         END DO
         ALL = .FALSE.
         GO TO 123
      ENDIF

      IF (NDEF /= 0) CALL SUMMRY (24)

      RETURN
      END SUBROUTINE GENRWF
