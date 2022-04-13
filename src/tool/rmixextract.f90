!
!***********************************************************************
!                                                                      *
PROGRAM extmix
!                                                                      *
!     Extract mixing coefficients and CSF's from files                 *
!     <name>.c, <name>.m / <name>.cm                                   *
!                                                                      *
!     Modified by G. Gaigalas                                   2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
   USE vast_kind_param, ONLY: DOUBLE
   USE eigv_C, ONLY: EAV, EVAL, EVEC
   USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   CHARACTER*100, line(3), g92mix*6, head*500
   CHARACTER*64 StrInput, basnam, from*1, suffix*3, filnam*69, dotc*2
   LOGICAL sort, getyn, first_of_the_block
   DATA dotc/'.c'/
   INTEGER :: nfmix, nfcsf, nfout, nfscratch
   DATA nfmix, nfcsf, nfout, nfscratch/20, 21, 22, 23/
   INTEGER :: ncfblk, nw, nvectot, ncftot, nvecsiz, nblock, layer
   INTEGER :: nevblk, nelec, nb, maxbub, icf, jcf, jblock, nmax
   INTEGER :: ivecdum, iaspa, iatjp, i, j, ip, IOS, IERR
   INTEGER :: icount, ans_all_states
   INTEGER, DIMENSION(:), pointer :: icount0
   REAL(DOUBLE) :: cutoff, bubmax
   INTEGER, DIMENSION(:, :), pointer :: iset0, jset0
   INTEGER, Allocatable, DIMENSION(:) :: iset, jset
!-----------------------------------------------
!
!***********************************************************************
!  Input conversation on the purpose, filenames, control parameters
!***********************************************************************

!     ...General description
   WRITE (*, *) '  RMIXEXTRACT'
   WRITE (*, *)
   WRITE (*, *) '  Extract and print mixing coefficients above a given'
   WRITE (*, *) '  cut-off. Resulting CSFs are written to screen and '
   WRITE (*, *) '  to rcsf.out.'
   WRITE (*, *)
   WRITE (*, *) '  Input files: name.c, name.(c)m'
   WRITE (*, *) '  Output file: rcsf.out'
   WRITE (*, *)

   !     ...Ask about which mode to run the code in

   WRITE (*, *) "  The extraction of CSF's can be done for each individual"
   WRITE (*, *) "  state separately, or across all states together."
   WRITE (*, *) "  The latter implies that a CSF will be extracted if has"
   WRITE (*, *) "  a mixing coefficient larger than the given cut-off in"
   WRITE (*, *) "  any of the eigenstates present in the mixing file."
   WRITE (*, *)
   WRITE (*, *) "- Do you want the extraction to be done"
   WRITE (*, *) "  for individual (0) or across all states (1)?: "
   READ (istdi, '(I5)') ans_all_states ! 1 = across all states, 0 = for individual states

   !     ...Asking the base name

123 WRITE (*, *) '- Give file name <name>.(c)m: '
   READ (istdi, '(A)') StrInput
   IF (LEN_TRIM(StrInput) .EQ. 0) THEN
      WRITE (istde, *) 'Blank line not acceptable, redo.'
      GOTO 123
   END IF

   basnam = ADJUSTL(StrInput)

!     ...Determing the suffix

234 WRITE (*, *) '- Are mixing coefficients from a CI calculation (y/n)? '
   READ (istdi, '(A)') StrInput
   StrInput = ADJUSTL(StrInput)
   from = StrInput(1:1)
   IF (from .EQ. 'y' .OR. from .EQ. 'Y') THEN
      suffix = '.cm'
   ELSEIF (from .EQ. 'n' .OR. from .EQ. 'N') THEN
      suffix = '.m'
   ELSE
      WRITE (*,*) StrInput(1:LEN_TRIM(StrInput)), &
         ' is not a valid choice, redo.'
      GOTO 234
   END IF

!     ...Form the mixing coefficient filename and check existence

   filnam = basnam(1:LEN_TRIM(basnam))//suffix
   INQUIRE (FILE=filnam, EXIST=sort) ! Borrow sort
   IF (.NOT. sort) THEN
      WRITE (*, *) 'File "', filnam(1:LEN_TRIM(filnam)), &
         '" does not exist, redo'
      GOTO 123
   END IF

!     ...The cut-off parameter

   WRITE (*, *) '- Enter the cut-off value for the coefficients ', &
      '[0--1]: '
   READ (istdi, *) cutoff
   cutoff = DABS(cutoff)

!     ...Sort or not

   PRINT *, '- Sort extracted CSFs by mixing coefficients (y/n)?'
   sort = getyn()

!***********************************************************************
!  Open mix file, check header
!***********************************************************************
   OPEN (nfmix, FILE=filnam, FORM='UNFORMATTED', STATUS='OLD' &
         , IOSTAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE (*, *) 'Failed to to open file: "', &
         filnam(1:LEN_TRIM(filnam)), '"'
      CLOSE (nfmix)
      STOP
   END IF

   READ (nfmix) g92mix
   IF (g92mix .NE. 'G92MIX') &
      STOP 'Not a mixing coefficient file!'

   READ (nfmix) nelec, ncftot, nw, nvectot, nvecsiz, nblock
   PRINT *
   PRINT '(2X,A,2X,I2,2X,A,2X,I8,2X,A,2X,I2,2X,A,2X,I2)', &
      'nblock = ', nblock, '  ncftot = ', ncftot, &
      '  nw = ', nw, '  nelec = ', nelec
   PRINT *

!***********************************************************************
!  Open CSF file
!***********************************************************************
   filnam = basnam(1:LEN_TRIM(basnam))//dotc
   OPEN (nfcsf, FILE=filnam, FORM='FORMATTED', STATUS='OLD' &
         , IOSTAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE (*, *) 'Failed to to open file: "', &
         filnam(1:LEN_TRIM(filnam)), '"'
      STOP
   END IF

   OPEN (nfout, FILE='rcsf.out', STATUS='UNKNOWN')

   ! The header lines of the CSF
   DO i = 1, 5
      READ (nfcsf, '(A)') head
      WRITE (nfout, '(A)') head(1:LEN_TRIM(head))
   END DO

!***********************************************************************
!  Loop over blocks. For each block
!    Read/write block info
!    Read CSF and write to a direct access file
!    Proceed if nevblk > 0
!       Allocate memory and read mix files
!***********************************************************************
   DO 432 jblock = 1, nblock

      first_of_the_block = .TRUE.

      READ (nfmix) nb, ncfblk, nevblk, iatjp, iaspa
      PRINT *, ('=', i=1, 75)
      PRINT '(2X,A,2X,I2,2X,A,2X,I8,3(2X,A,2X,I2))', &
         'nb = ', nb, 'ncfblk = ', ncfblk, &
         'nevblk = ', nevblk, '2J+1 = ', iatjp, 'parity = ', iaspa
      WRITE (istde, '(2X,A,2X,I2,2X,A,2X,I8,3(2X,A,2X,I2))') &
         'nb = ', nb, 'ncfblk = ', ncfblk, &
         'nevblk = ', nevblk, '2J+1 = ', iatjp, 'parity = ', iaspa
      PRINT *, ('=', i=1, 75)
      IF (jblock .NE. nb) STOP 'jblock .NE. nb'

      CALL iocsf(nfcsf, nfscratch, jblock, ncfblk, line)

      IF (nevblk .LE. 0) GOTO 432

      Allocate (eval(1:nevblk), stat=ierr)
      IF (ierr /= 0) STOP " Not enough memory for eval!"
      Allocate (evec(1:nevblk*ncfblk), stat=ierr)
      IF (ierr /= 0) STOP " Not enough memory for evec!"
      IF (ans_all_states == 0) THEN
         Allocate (icount0(1:nevblk), stat=ierr)
         IF (ierr /= 0) STOP " Not enough memory for icount0!"
         Allocate (iset0(1:ncfblk, 1:nevblk), stat=ierr)
         IF (ierr /= 0) STOP " Not enough memory for iset0!"
      ELSE
         Allocate (iset(1:ncfblk), stat=ierr)
         IF (ierr /= 0) STOP " Not enough memory for iset"
      END IF

      READ (nfmix) (ivecdum, i=1, nevblk)
      READ (nfmix) eav, (eval(i), i=1, nevblk)
      READ (nfmix) ((evec(i+(j-1)*ncfblk), i=1, ncfblk), j=1, nevblk)

      IF (ans_all_states == 0) THEN
!        ...Find the "OR" set
         icount0 = 0
         DO icf = 1, ncfblk
            DO ip = 1, nevblk
               IF (DABS(evec(icf + (ip - 1)*ncfblk)) > cutoff) THEN
                  icount0(ip) = icount0(ip) + 1
                  iset0(icount0(ip), ip) = icf
               END IF
            END DO
         END DO

         DO ip = 1, nevblk
            if (ip == 1) then
               PRINT *,                                                &
               'Average Energy = ',eav,'    ncf_reduced = ',icount0(ip)
            else
               PRINT *, '                    Eigenvector =', ip,       &
               '  ncf_reduced = ', icount0(ip)
            end if
            nmax = max(nmax, icount0(ip))
         END DO

!        ...Make a copy of the original set which is to be altered if
!           sorted
         IF (sort) THEN
            allocate (jset0(1:nmax, 1:nevblk))
            DO j = 1, nevblk
               DO i = 1, icount0(j)
                  jset0(i, j) = iset0(i, j)
               END DO
            END DO
         END IF

         DO 1 ip = 1, nevblk
            layer = (ip - 1)*ncfblk
            eval(ip) = eval(ip) + eav

            IF (sort) THEN
               DO i = 1, icount0(ip)
                  icf = iset0(i, ip)
                  maxbub = i
                  bubmax = ABS(evec(icf + layer))
                  DO j = i + 1, icount0(ip)
                     jcf = iset0(j, ip)
                     IF (ABS(evec(jcf + layer)) .GT. bubmax) THEN
                        maxbub = j
                        bubmax = ABS(evec(jcf + layer))
                     END IF
                  END DO

                  iset0(i, ip) = iset0(maxbub, ip)
                  iset0(maxbub, ip) = icf
               END DO
            END IF

            PRINT *
            PRINT *, 'Energy = ', eval(ip), '    Coefficients and CSF :'
            PRINT *
            DO i = 1, icount0(ip)
               icf = iset0(i, ip)
               write (*, '(i12,f11.6)') i, evec(icf + layer)
               READ (nfscratch, REC=icf) line
               PRINT *, line(1) (1:LEN_TRIM(line(1)))
               PRINT *, line(2) (1:LEN_TRIM(line(2)))
               PRINT *, line(3) (1:LEN_TRIM(line(3)))

               IF (first_of_the_block) THEN
                  WRITE (nfout, '(A)') line(1) (1:LEN_TRIM(line(1)))
                  WRITE (nfout, '(A)') line(2) (1:LEN_TRIM(line(2)))
                  WRITE (nfout, '(A)') line(3) (1:LEN_TRIM(line(3)))
               END IF

!              ...Recover original set
               IF (sort) iset0(i, ip) = jset0(i, ip)
            END DO

         first_of_the_block = .FALSE.

   1     CONTINUE
         deallocate (iset0)
            IF (sort) deallocate (jset0)
      ELSE

!        ...Find the "OR" set
         icount = 0
         DO icf = 1, ncfblk
            DO ip = 1, nevblk
               IF (ABS(evec(icf + (ip - 1)*ncfblk)) .GT. cutoff) THEN
                  icount = icount + 1
                  iset(icount) = icf
                  EXIT
               END IF
            END DO
         END DO

!        ...Make a copy of the original set which is to be altered if
!           sorted
         IF (sort) THEN
            allocate (jset(1:icount))
            DO i = 1, icount
               jset(i) = iset(i)
            END DO
         END IF

         PRINT *, 'Average Energy = ', eav, '    ncf_reduced = ', icount

         DO 2 ip = 1, nevblk
            layer = (ip - 1)*ncfblk
            eval(ip) = eval(ip) + eav

            IF (sort) THEN
               DO i = 1, icount
                  icf = iset(i)
                  maxbub = i
                  bubmax = ABS(evec(icf + layer))
                  DO j = i + 1, icount
                     jcf = iset(j)
                     IF (ABS(evec(jcf + layer)) .GT. bubmax) THEN
                        maxbub = j
                        bubmax = ABS(evec(jcf + layer))
                     END IF
                  END DO

                  iset(i) = iset(maxbub)
                  iset(maxbub) = icf
               END DO
            END IF

            PRINT *
            PRINT *, 'Energy = ', eval(ip), '    Coefficients and CSF :'
            PRINT *

            DO i = 1, icount
               icf = iset(i)
               WRITE (*, '(i12,f11.6)') i, evec(icf + layer)
               READ (nfscratch, REC=icf) line
               PRINT *, line(1) (1:LEN_TRIM(line(1)))
               PRINT *, line(2) (1:LEN_TRIM(line(2)))
               PRINT *, line(3) (1:LEN_TRIM(line(3)))

               IF (first_of_the_block) THEN
                  WRITE (nfout, '(A)') line(1) (1:LEN_TRIM(line(1)))
                  WRITE (nfout, '(A)') line(2) (1:LEN_TRIM(line(2)))
                  WRITE (nfout, '(A)') line(3) (1:LEN_TRIM(line(3)))
               END IF

!              ...Recover original set
               IF (sort) iset(i) = jset(i)
            END DO

            first_of_the_block = .FALSE.

   2     CONTINUE
         deallocate (iset)
         IF (sort) deallocate (jset)
      END IF

      IF (jblock .LT. nblock) WRITE (nfout, '(A)') ' *'

      deallocate (eval)
      deallocate (evec)

      CLOSE (nfscratch, STATUS='DELETE')

432 CONTINUE

   CLOSE (nfmix)
   CLOSE (nfcsf)
   CLOSE (nfout)

   WRITE (*, *) 'RMIXEXTRACT: Execution complete.'

CONTAINS

   SUBROUTINE iocsf(nfcsf, nfscratch, jblock, ncfblk, line)
      IMPLICIT NONE
      INTEGER nfcsf, nfscratch, jblock, ncfblk, icf, i
      CHARACTER*(*) line(3), star*2

      OPEN (nfscratch, STATUS='SCRATCH', ACCESS='DIRECT', &
            RECL=300)

      DO icf = 1, ncfblk
         READ (nfcsf, '(A)') line(1)
         READ (nfcsf, '(A)') line(2)
         READ (nfcsf, '(A)') line(3)
         WRITE (nfscratch, REC=icf) line
      END DO

      READ (nfcsf, '(A)', END=123) star
      IF (star .NE. ' *') STOP ' Error while reading CSFs!'

123   CONTINUE

      RETURN
   END SUBROUTINE iocsf

END PROGRAM extmix
