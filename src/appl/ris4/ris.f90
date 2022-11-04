!***********************************************************************
!                                                                      *
      PROGRAM RIS
!                                                                      *
!   Entry routine for RIS. Controls the entire computation.       *
!                                                                      *
!   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
!               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM                *
!                        STRSUM                                        *
!               [NJGRAF]: FACTT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: 17 Jan 1996   *
!   Modify  by Gediminas Gaigalas                        26 Oct 2009   *
!   Modified by J. Ekman                                 25 Nov 2013   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE default_C
      USE iounit_C
      USE debug_C,         ONLY: CUTOFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcsla_I
      USE getmixblock_I
      USE strsum_I
      USE factt_I
      USE ris_cal_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: DR2
      LOGICAL   :: YES
      CHARACTER :: NAME*24
      INTEGER :: K, NCI, ncore_not_used, NOPAR
!-----------------------------------------------
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
      CUTOFF = 1.0D-10
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'RIS: Execution begins ...'
!      CALL STARTTIME (ncount1, 'RIS')

      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'Default settings?'
      YES = GETYN ()
      WRITE (ISTDE, *)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

    9 WRITE (ISTDE, *) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE (ISTDE, *) 'Names may not start with a blank'
         GOTO 9
      ENDIF
      PRINT *

      WRITE (ISTDE, *) 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
      WRITE (ISTDE, *)

!
!   Check compatibility of plant substitutions
!
!GG      CALL CHKPLT
!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      CALL SETDBG
!
!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Open the  .i  file
!
      CALL SETSUM(NAME,NCI)
!
!   Open, check, load data from, and close, the  .csl  file
!
      CALL SETCSLA(NAME,ncore_not_used)
!
!   Get the remaining information
!
      CALL GETSMD(NAME)
!
!   Get the eigenvectors
!
!      PRINT *, 'Block format?'
!      YES = GETYN ()
!      PRINT *

!      IF (YES) THEN
         CALL GETMIXBLOCK(NAME,NCI)
!      ELSE
!         IF (NCI.EQ.0) THEN
!            CALL GETMIXC(NAME)
!         ELSE
!            CALL GETMIXA(NAME)
!         ENDIF
!      ENDIF
!
!   Append a summary of the inputs to the  .sum  file
!
!      CALL STRSUM
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the RIS calculation
!
      CALL RIS_CAL(NAME)
!
!   Print completion message
!
      WRITE (ISTDE, *)
!      CALL STOPTIME (ncount1, 'RIS')
      WRITE (ISTDE, *) 'RIS: Execution complete.'
!
      STOP
      END  PROGRAM RIS
