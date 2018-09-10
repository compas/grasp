!***********************************************************************
!                                                                      *
      SUBROUTINE GETSCD(EOL, IDBLK, ISOFILE, RWFFILE) 
!                                                                      *
!   Interactively determines the data governing the SCF problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETRSL, GETYN, NUCPOT, RADGRD,        *
!                        SETISO, SETQIC, SETRWF.                       *
!               [RSCF92]: GETALD, GETOLD
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Dec 1992   *
!
! Xinghong He 98-08-06
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE 
      USE parameter_def,    ONLY: NNNW, NNNP
      USE COUN_C 
      USE damp_C, ONLY: odamp, cdamp
      USE DEF_C 
      USE default_C
      USE fixd_C, ONLY: lfix
      USE iounit_C
      USE grid_C, ONLY: h, hp, n, rnt
      USE invt_C, ONLY: noinvt
      USE mcpb_C
      USE mpi_C
      USE node_C, ONLY: nnodep
      USE npar_C
      !USE npot_C
      USE ORB_C , ONLY: nak, nh, np, nw, PED
      USE orthct_C
      USE scf_C, ONLY: scnsty, method
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE setiso_I 
      USE setqic_I 
      USE radgrd_I 
      USE nucpot_I 
      USE setrwfa_I 
      USE getald_I 
      USE getold_I 
      USE convrt_I 
      USE getrsl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(OUT) :: EOL 
      CHARACTER  :: ISOFILE*(*) 
      CHARACTER  :: RWFFILE*(*) 
      CHARACTER  :: IDBLK(*)*8 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*6, PARAMETER :: MYNAME = 'GETSCD' 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQADUM 
      INTEGER, DIMENSION(NNNW) :: INDX 
      INTEGER :: I, J, IEND, IBEG, LENTH, NSUBS, LOC 
      REAL(DOUBLE) :: ODAMPU, CDAMPU 
      LOGICAL :: YES 
      CHARACTER :: RECORD*80, CNUM*20 
!-----------------------------------------------
!
!   Open, check, load data from, and close the  .iso  file
!
      CALL SETISO (ISOFILE)
!
!   Set default speed of light and grid parameters
!
      C = CVAC 
      IF (NPARM == 0) THEN 
         RNT = EXP((-65.D0/16.D0))/Z 
         H = 0.5D0**4 
         N = MIN(220,NNNP) 
      ELSE 
!         default comes here
!CFF     .. should be Z-dependent
         RNT = 2.0D-06/Z 
         H = 5.D-2 
         N = NNNP 
      ENDIF 
      HP = 0.D0 
!
      IF (NDEF /= 0) THEN 
 
         WRITE (ISTDE,'(A)')'Change the default speed of',&
            ' light or radial grid parameters?  (y/n) ' 
         YES = GETYN() 
         IF (YES) THEN 
            WRITE (istde,*) 'Speed of light = ',CVAC,';', ' revise ?'
            YES = GETYN() 
            IF (YES) THEN 
               WRITE (ISTDE, *) 'Enter the revised value:' 
               READ (5, *) C 
            ENDIF 
!
!   Determine the parameters controlling the radial grid
!
            WRITE (ISTDE, *) 'The default radial grid parameters for ',&
               'this case are:' 
            WRITE (ISTDE, *) ' RNT = ', RNT 
            WRITE (ISTDE, *) ' H = ', H 
            WRITE (ISTDE, *) ' HP = ', HP 
            WRITE (ISTDE, *) ' N = ', N 
            WRITE (ISTDE, *) ' revise these values?' 
            YES = GETYN() 
            IF (YES) THEN 
               WRITE (ISTDE, *) 'Enter RNT:' 
               READ (5, *) RNT 
               WRITE (ISTDE, *) 'Enter H:' 
               READ (5, *) H 
               WRITE (ISTDE, *) 'Enter HP:' 
               READ (5, *) HP 
               WRITE (ISTDE, *) 'Enter N:' 
               READ (5, *) N 
!b Revised grid
               WRITE (istde,*) 'Revised RNT = ', RNT
               WRITE (istde,*) 'Revised H   = ', H
               WRITE (istde,*) 'Revised HP  = ', HP
               WRITE (istde,*) 'Revised N   = ', N
            ENDIF 
         ENDIF 
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      ACCY = H**6
!
!b
!b read ACCY on input
!b
         WRITE (istde,*) 'Revise the default ACCY = ', ACCY
         YES = GETYN ()
            IF (YES) THEN
               WRITE (istde,*) 'Enter ACCY:'
               READ *, ACCY
!b Revised ACCY
               WRITE (istde,*) 'Revised ACCY = ', ACCY
         ENDIF
      ENDIF
!
!   Set up the coefficients for the numerical procedures
!
      CALL SETQIC 
!
!   Generate the radial grid and all associated arrays
!
      CALL RADGRD 
!
!   Generate $- r \times V_ (r)$
!
      CALL NUCPOT 
!
!   Load the subshell radial wavefunction estimates
!
      CALL SETRWFA (RWFFILE) 
!
!   Set some defaults
!
      THRESH = 0.05D0 
 
!      IORDER(I) = I    ! Completely determined in GETOLD
      METHOD(:NW) = 3 
      NOINVT(:NW) = .TRUE. 
!CFF      ODAMP(:NW) = 1.D0 
      ODAMP(:NW) = 0.0D0 
      PED(:NW)   = 0.0D0
      SCNSTY(:NW) = 0.0D0 
 
      WHERE (NAK(:NW) < 0)  
         NNODEP(:NW) = NP(:NW) + NAK(:NW) 
      ELSEWHERE 
         NNODEP(:NW) = NP(:NW) - NAK(:NW) - 1 
      END WHERE 
 
      IF (DIAG) THEN 
         EOL = .FALSE. 
         CALL GETALD          ! EOL type calculation, 
                              ! H(DC) will not be diagonalised
      ELSE IF (LFORDR) THEN 
         EOL = .TRUE. 
         CALL GETOLD (IDBLK)  ! EOL type calculation 
      ELSE 
!GG         WRITE (ISTDE, '(A)') 'EOL type calculation?  (y/n) ' 
         EOL = .true.
!GG         EOL = GETYN() 
         IF (EOL) THEN 
            CALL GETOLD (IDBLK) 
         ELSE 
            CALL GETALD       ! EAL type calculation, 
                              ! H(DC) will not be diagonalised
         ENDIF 
      ENDIF 
 
      WRITE (ISTDE, *) 'Enter the maximum number of SCF cycles:' 
      READ (*, *) NSCF 

      IF (NDEF.EQ.0) THEN
         WRITE(734,*) NSCF,'! Number of SCF cycles'
      END IF
!
!   Allow the user to modify other defaults
!
      IF (NDEF /= 0) THEN 
         WRITE (ISTDE, '(A)') 'Modify other defaults?  (y/n) ' 
         YES = GETYN() 
      ELSE 
         YES = .FALSE. 
      ENDIF 
 
      IF (.NOT.YES) RETURN  
!=======================================================================
! From here to end, "other defaults" are handled. For simplicity
! We'll let node-0 do the job and then broadcast results to all
! nodes.
!=======================================================================
                  !-------------------------------------------
      IF (MYID == 0) THEN                        ! This is a _big_ IF 
                  !-------------------------------------------
!
!   THRESH
!
         WRITE (ISTDE, *)'An oscillation in the large-component of the ',&
            'radial wavefunction is diregarded' 
         WRITE (ISTDE, *) 'for the purposes of node counting if its ', &
            'amplitude is less than 1/20 the' 
         WRITE (ISTDE, *) 'maximum amplitude.   Revise this?' 
         YES = GETYN() 
         IF (YES) THEN 
    3       CONTINUE 
            WRITE (ISTDE, *) 'Enter the new threshold value:' 
            READ (*, *) THRESH 
            IF (THRESH <= 0.D0) THEN 
               WRITE (ISTDE, *) MYNAME, ': This must exceed 0;' 
               GO TO 3 
            ENDIF 
         ENDIF 
         YES = .FALSE.
!
!   METHOD
!
 
! Piece only for printing...
 
         DO I = 1, 4 
 
            DO J = 1, NW 
               IF (.NOT.(METHOD(J)==I .AND. .NOT.LFIX(J))) CYCLE  
               WRITE (ISTDE, *) 'Method ', I, ' is used for ', &
                  'integrating the radial differential ', &
                  'equation for subshells' 
               GO TO 9 
            END DO 
 
            CYCLE  
 
    9       CONTINUE 
            IEND = 0 
            DO J = 1, NW 
               IF (METHOD(J)==I .AND. .NOT.LFIX(J)) THEN 
                  IBEG = IEND + 1 
                  IEND = IBEG 
                  RECORD(IBEG:IEND) = ' ' 
                  CALL CONVRT (NP(J), CNUM, LENTH) 
                  IBEG = IEND + 1 
                  IEND = IBEG + LENTH - 1 
                  RECORD(IBEG:IEND) = CNUM(1:LENTH) 
                  IBEG = IEND + 1 
                  IF (NAK(J) < 0) THEN 
                     IEND = IBEG 
                     RECORD(IBEG:IEND) = NH(J)(1:1) 
                  ELSE 
                     IEND = IBEG + 1 
                     RECORD(IBEG:IEND) = NH(J)(1:2) 
                  ENDIF 
               ENDIF 
               IF (IEND <= 76) CYCLE  
               WRITE (ISTDE, *) RECORD(1:IEND) 
               IEND = 0 
            END DO 
            IF (IEND<=0 .OR. MYID/=0) CYCLE  
            WRITE (ISTDE, *) RECORD(1:IEND) 
         END DO 
 
! Reads user inputs and fills array indx(1:nsubs) where nsubs itself
! is an output from getrsl.
! METHOD(1:4) is the only output. indx() and nsubs are discarded
 
         WRITE (ISTDE, *) 'Select a different integration method for ', &
            'any subshell radial wavefunction?' 
         YES = GETYN() 
         IF (YES) THEN 
            DO I = 1, 4 
               WRITE (ISTDE, *) 'Method ', I, ':' 
               CALL GETRSL (INDX, NSUBS) 
               DO J = 1, NSUBS 
                  LOC = INDX(J) 
                  IF (LFIX(LOC)) CYCLE  
                  METHOD(LOC) = I 
               END DO 
            END DO 
         ENDIF 
!
!   NOINVT
!
         WRITE (ISTDE, *) 'The first oscillation of the large component' 
 
         DO I = 1, NW 
            IF (.NOT.(NOINVT(I) .AND. .NOT.LFIX(I))) CYCLE  
            WRITE (ISTDE, *) 'of the following radial wavefunctions ', &
               'will be required to be positive' 
            GO TO 15 
         END DO 
 
         WRITE (ISTDE, *) 'of all radial wavefunctions will be required ', &
            'to be positive.   Revise this?' 
         YES = GETYN() 
         GO TO 17 
   15    CONTINUE 
         IEND = 0 
         DO I = 1, NW 
            IF (NOINVT(I) .AND. .NOT.LFIX(I)) THEN 
               IBEG = IEND + 1 
               IEND = IBEG 
               RECORD(IBEG:IEND) = ' ' 
               CALL CONVRT (NP(I), CNUM, LENTH) 
               IBEG = IEND + 1 
               IEND = IBEG + LENTH - 1 
               RECORD(IBEG:IEND) = CNUM(1:LENTH) 
               IBEG = IEND + 1 
               IF (NAK(I) < 0) THEN 
                  IEND = IBEG 
                  RECORD(IBEG:IEND) = NH(I)(1:1) 
               ELSE 
                  IEND = IBEG + 1 
                  RECORD(IBEG:IEND) = NH(I)(1:2) 
               ENDIF 
            ENDIF 
            IF (IEND <= 76) CYCLE  
            WRITE (ISTDE, *) RECORD(1:IEND) 
            IEND = 0 
         END DO 
         IF (IEND > 0) WRITE (ISTDE, *) RECORD(1:IEND) 
         WRITE (ISTDE, *) 'Revise this?' 
         YES = GETYN() 
   17    CONTINUE 
         IF (YES) THEN 
            WRITE (ISTDE, *) 'Suppressing enforcement of positive first ', &
               'oscillation:' 
            CALL GETRSL (INDX, NSUBS) 
            DO I = 1, NSUBS 
               LOC = INDX(I) 
               IF (LFIX(LOC)) CYCLE  
               NOINVT(LOC) = .TRUE. 
            END DO 
         ENDIF 
!
!   ODAMP
!
         DO I = 1, NW 
            IF (.NOT.(ODAMP(I)/=0.D0 .AND. .NOT.LFIX(I))) CYCLE  
            WRITE (ISTDE, *) 'Subshell accelerating parameters have ', &
               'been set.   Revise these?' 
            YES = GETYN() 
            GO TO 20 
         END DO 
         WRITE (ISTDE, *) 'Set accelerating parameters for subshell ', &
            'radial wavefunctions?' 
         YES = GETYN() 
   20    CONTINUE 
         IF (YES) THEN 
            WRITE (ISTDE, *) 'Different accelerating parameters for ', &
               'different subshell radial wavefunction?' 
            YES = GETYN() 
            IF (YES) THEN 
   21          CONTINUE 
               WRITE (ISTDE, *) 'Enter an accelerating parameter' 
               WRITE (ISTDE, *) ' (0< ODAMP < 1 allows ODAMP to be ', &
                  'reduced as convergence is approached;' 
               WRITE (ISTDE, *) ' -1 < ODAMP < 0 implies |ODAMP| is ', &
                  'held constant):' 
               READ (*, *) ODAMPU 
               IF (ABS(ODAMPU)==0.D0 .OR. ABS(ODAMPU)>=1.D0) THEN 
                  WRITE (ISTDE, *) MYNAME, ': Value out of range ...' 
                  GO TO 21 
               ELSE 
                  CALL GETRSL (INDX, NSUBS) 
                  DO I = 1, NSUBS 
                     LOC = INDX(I) 
                     IF (LFIX(LOC)) CYCLE  
                     ODAMP(LOC) = ODAMPU 
                  END DO 
               ENDIF 
            ELSE 
   23          CONTINUE 
               WRITE (ISTDE, *) 'Enter the accelerating parameter' 
               WRITE (ISTDE, *) ' (0< ODAMP < 1 allows ODAMP to be ', &
                  'reduced as convergence is approached;' 
               WRITE (ISTDE, *) ' -1 < ODAMP < 0 implies |ODAMP| is ', &
                  'held constant):' 
               READ (*, *) ODAMPU 
               IF (ABS(ODAMPU)==0.D0 .OR. ABS(ODAMPU)>=1.D0) THEN 
                  WRITE (ISTDE, *) MYNAME, ': Value out of range ...' 
                  GO TO 23 
               ELSE 
                  WHERE (.NOT.LFIX(:NW))  
                     ODAMP(:NW) = ODAMPU 
                  END WHERE 
               ENDIF 
            ENDIF 
         ENDIF 
!
!   CDAMP
!
         WRITE (ISTDE, *) 'Set accelerating parameters for the ', &
            'eigenvectors?' 
         YES = GETYN() 
         IF (YES) THEN 
            WRITE (ISTDE, *) 'Different accelerating parameters ', &
               'for each eigenvector?' 
            YES = GETYN() 
            IF (YES) THEN 
               WRITE (ISTDE, *) 'Enter an accelerating parameter for' 
               CALL CONVRT (NCMIN, RECORD, LENTH) 
               WRITE (ISTDE, *) ' each of the '//RECORD(1:LENTH)//' levels :' 
               READ (*, *) (CDAMP(I),I=1,NCMIN) 
            ELSE 
               WRITE (ISTDE, *) 'Enter the accelerating parameter:' 
               READ (*, *) CDAMPU 
               CDAMP(:NCMIN) = CDAMPU 
            ENDIF 
         ENDIF 
!
!   NSIC
!
         WRITE (ISTDE, *) 'Following the improvement of each of the ', &
            'subshell radial wavefunctions in turn, ' 
         WRITE (ISTDE, *) 'the ', NSIC, ' least self-consistent', &
            ' functions will be improved at the' 
         WRITE (ISTDE, *) 'end of the first SCF cycle. Revise this ', &
            'setting?' 
         YES = GETYN() 
         IF (YES) THEN 
            WRITE (ISTDE, *) 'Enter the number of additional ', 'improvements:' 
            READ (*, *) NSIC 
         ENDIF 
!
!   NSOLV
!
         WRITE (ISTDE, *) 'The maximum number of cycles in attempting ', &
            'to solve each radial equation is ' 
         WRITE (ISTDE, *) NSOLV, ' times the principal quantum', &
            ' number of the radial' 
         WRITE (ISTDE, *) 'wave-function to be estimated.   ', &
            'Revise this setting?' 
         YES = GETYN() 
         IF (YES) THEN 
            WRITE (ISTDE, *) 'Enter the factor that multiplies the ', &
               'principal quantum number:' 
            READ (*, *) NSOLV 
         ENDIF 
!
!   Orthogonalisation
!
         IF (ORTHST) THEN 
            WRITE (ISTDE, *) 'Subshell radial wavefunctions will be ', &
               'Schmidt orthogonalised immediately' 
            WRITE (ISTDE, *) 'following their estimation to all ', &
               'functions with poorer self-consistency.' 
            WRITE (ISTDE, *) ' Revise this?' 
            YES = GETYN() 
            IF (YES) ORTHST = .FALSE. 
         ELSE 
            WRITE (ISTDE, *) 'Subshell radial wavefunctions will be ', &
               'Schmidt orthogonalised at the end of' 
            WRITE (ISTDE, *) 'each SCF cycle.   Revise this?' 
            YES = GETYN() 
            IF (YES) ORTHST = .TRUE. 
         ENDIF 
                  !-------------------------------------------
      ENDIF                                      ! end of the _big_ IF 
                  !-------------------------------------------
      RETURN  
      END SUBROUTINE GETSCD 
