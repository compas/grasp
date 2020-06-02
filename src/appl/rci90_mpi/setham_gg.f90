!***********************************************************************
!                                                                      *
      SUBROUTINE SETHAM (myid, nprocs, jblock, ELSTO,ICSTRT, nelmntt, &
                         atwinv,slf_en)
!                                                                      *
!   Sets up the Hamiltonian matrix and determines the average energy.  *
!
!   Serial I/O moved out; able to run on single/multiple processors
!   For this purpose a new common /setham_to_genmat2/ is created
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
!               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: NNNP, KEYORB
      USE memory_man
!cjb iso_fortran_env :: output_unit, error_unit
      use, intrinsic :: iso_fortran_env
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE bcore_C, ONLY: icore
      USE bilst_C
      USE buffer_C
      USE coeils_C
      USE decide_C
      USE debug_C
      USE def_C
      USE eigv_C
      USE iccu_C
      USE grid_C
      USE hmat_C
      USE keilst_C
      USE ncdist_C
      USE orb_C, ONLY: ncf,  nw
      USE prnt_C
      USE stat_C
      USE stor_C
      USE tatb_C
      USE vinlst_C
      USE vpilst_C
      USE wave_C
      USE cteilsrk_C
      USE blim_c
      USE where_C
      USE eigvec1_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE iabint_I
      USE brint1_I
      USE brint2_I
      USE brint3_I
      USE brint4_I
      USE brint5_I
      USE brint6_I
      USE convrt_I
      USE ichop_I
      USE keint_I
      USE rkintc_I
      USE vint_I
      USE vpint_I
      IMPLICIT NONE
      EXTERNAL BREID,CORD
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER(LONG)       :: NELMNTT
      INTEGER             :: JBLOCK, ICSTRT
      INTEGER, INTENT(IN) :: MYID, NPROCS
      REAL(DOUBLE)        :: ELSTO, ATWINV
      REAL(DOUBLE), DIMENSION(*) :: slf_en
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      REAL(DOUBLE), DIMENSION(NNNW) :: tshell
      REAL(DOUBLE) :: tgrl1, tgrl2, tegral

      INTEGER, PARAMETER :: KEY = KEYORB
!
!     Matrix elements smaller than CUTOFF are not accumulated
!
!cjb cutoff is use associated and cannot be redeclared
!     REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!cjb CUTOFF = 1.0D-12 below
!cjb
      INTEGER :: ipi, ipj, inc1, inc2, kt, ipt, incor, ncoec, nctec, &
                 i, j, nmcbp, ncore, ic, nelc, irstart, ir, ia, ib,  &
                 itype, nctei, iia
       REAL(DOUBLE) :: elemnt, precoeff, tcoeff, vcoeff, contr
!-----------------------------------------------------------------------
      nelmnt = nelmntt
!cjb
!cjb  CUTOFF = 1.0D-12
      CUTOFF = 1.0D-12

      ATWINV = 1.D0/EMN

      IF (IPRERUN .EQ. 2) THEN
         DO IPI = 1,NVEC
            DO IPJ = 1,NCF
               WRITE (*,*) IPI,IPJ,EVEC1(IPJ+(IPI-1)*NCF)
            ENDDO
         ENDDO
      ENDIF
!
!   Allocate storage to arrays in COMMON/BUFFER/; these are
!   used for the Coulomb and transverse two-electron integrals
!
      CALL ALCBUF (1)

!     ...Locals
      CALL alloc (emt, ncf, 'EMT','SETHAM' )
      CALL alloc (irow, ncf, 'IROW', 'SETHAM')
!
      INC1 = 1
      INC2 = 1
!
!   Initialisations for contributions from the Dirac-Coulomb
!   operator
!
      KT  = 0
      IPT = 1
!
      INCOR = 1

      NCOEC = 0
!
      NCTEC   = 0

      IF (LTRANS) THEN

!        ...Initialisations for transverse interaction correction
         DO 2 I = 1, NW
            ICORE(I) = 0
            DO J = 1, NCF
               IF (ICHOP (I,J) .LE. 0) GOTO 2
            ENDDO
            ICORE(I) = 1
    2    CONTINUE

         NMCBP = 0
         NCORE = 0
      ENDIF

! Loop over rows of the Hamiltonian matrix - distributed

      DO 10 ic = icstrt, ncf, nprocs

         NELC = 0    ! counter - Number of non-zeros of this row

! Loop over columns of the current row

         irstart = 1
         DO 85 IR = irstart, IC

! PER
            IF (LFORDR .AND. (IR .GT. ICCUT(1))) THEN
               IF (IR.NE.IC) CYCLE
            END IF
! PER

            ELEMNT = 0.D0     ! accumulates various contributions to H

!
!   Generate the integral list for the matrix element of the
!   one-body operators
!
            IF (IPRERUN .EQ. 1) THEN
               INC1 = 0
               INC2 = 0
               IF (IC.LE.NCSFPRE .OR. IC.EQ.IR) THEN
                  INC1 = 1
               ENDIF
            ENDIF

            IF (IPRERUN .EQ. 2) THEN
!
!   Diagonal elements are always included
!   Off diagonal elements are included only if the
!   products of the weights from the prerun are larger
!   than the cutoff.
!
               IF (IC .EQ. IR) THEN
                  INC1 = 1
                  INC2 = 1
               ELSE
                  INC1 = 0
                  INC2 = 0
               ENDIF
               DO IPI = 1,NVEC
                  PRECOEFF =                                           &
                   DABS(EVEC1(IC+(IPI-1)*NCF)*EVEC1(IR+(IPI-1)*NCF))
                  IF (PRECOEFF .GT. COEFFCUT1) INC1 = 1
                  IF (PRECOEFF .GT. COEFFCUT2) INC2 = 1
               ENDDO
            ENDIF

!            ...INC1.EQ.1 ------------>
         IF (INC1 .EQ. 1) THEN   !inc1 is always 1 without PRE-RUN
           CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
!
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
         IF (IA .NE. 0) THEN
            IF (IA .EQ. IB) THEN
               DO IA = 1,NW
                  iia = ia
                  TCOEFF = TSHELL(IA)
                  IF (ABS (TCOEFF) .GT. CUTOFF) THEN
                     NCOEC = NCOEC + 1
                     CALL IABINT (IIA, IIA, TEGRAL)
                        !------------------------
                     ELEMNT = ELEMNT + TEGRAL*TCOEFF
                     IF (LNMS) THEN
                        CALL KEINT (IIA,IIA,TEGRAL)
                        !------------------------
                        ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
                     ENDIF
                     IF (LVP) THEN
                        CALL VPINT (IIA, IIA, TEGRAL)
                        !------------------------
                        ELEMNT = ELEMNT + TEGRAL*TCOEFF
                     ENDIF
                  ENDIF
               ENDDO
            ELSE
               TCOEFF = DBLE(TSHELL(1))
               IF (DABS (TCOEFF) .GT. CUTOFF) THEN
                  NCOEC = NCOEC + 1
                  CALL IABINT (IA, IB, TEGRAL)
                        !------------------------
                  ELEMNT = ELEMNT + TEGRAL*TCOEFF
                  IF (LNMS) THEN
                     CALL KEINT (IA, IB, TEGRAL)
                        !------------------------
                     ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
                  ENDIF
                  IF (LVP) THEN
                     CALL VPINT (IA, IB, TEGRAL)
                        !------------------------
                     ELEMNT = ELEMNT + TEGRAL*TCOEFF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
         IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
         NVCOEF = 0
!
         CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
!
         DO 7 I = 1, NVCOEF
            VCOEFF = DBLE(COEFF(I))
            IF (DABS (VCOEFF) .GT. CUTOFF) THEN
               NCTEC = NCTEC + 1
               IF (LSMS) THEN
                  IF (LABEL(5,I) .EQ. 1) THEN
                     CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                     CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                     ELEMNT = ELEMNT - TGRL1*TGRL2*ATWINV*VCOEFF
                  ENDIF
               ENDIF
               CALL RKINTC (LABEL(1,I), LABEL(2,I),                  &
                            LABEL(3,I), LABEL(4,I),                  &
                            LABEL(5,I), TEGRAL)
               ELEMNT = ELEMNT + TEGRAL*VCOEFF
            ENDIF
    7    CONTINUE
!
         IBUG1 = 0

         ENDIF  !inc1 is always 1 without PRE-RUN
!            ...INC1.EQ.1 <------------
!***********************************************************************
!            ...LTRANS .AND. (INC2.EQ.1) ------------>
         IF (LTRANS .AND. (INC2.EQ.1)) THEN
            !IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
!
            NVCOEF = 0
!
            CALL RKCO_GG (IC, IR, BREID, 1, 2)
!
            DO 8 I = 1, NVCOEF
               IF (DABS (COEFF(I)) > CUTOFF) THEN
                  NMCBP = NMCBP + 1
                  ITYPE = ABS (LABEL(6,I))
                  IF (ITYPE == 1) THEN
                     CALL BRINT1 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE == 2) THEN
                     CALL BRINT2 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE == 3) THEN
                     CALL BRINT3 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE == 4) THEN
                     CALL BRINT4 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE == 5) THEN
                     CALL BRINT5 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE == 6) THEN
                     CALL BRINT6 (LABEL(1,I), LABEL(2,I),               &
                                  LABEL(3,I), LABEL(4,I),               &
                                  LABEL(5,I), TEGRAL)
                  ENDIF
                  CONTR = COEFF(I)*TEGRAL
                  IF (LABEL(6,I) > 0) THEN
                     ELEMNT = ELEMNT + CONTR
                  ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     NCORE = NCORE + 1
                     ELSTO = ELSTO + CONTR
                  ENDIF
               ENDIF
    8       CONTINUE
!
            IBUG1 = 0
!
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
!            IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO
!
            !ENDIF   !inc2 is always 1 without PRE-RUN
         ENDIF
!            ...LTRANS .AND. (INC2.EQ.1) <------------
!***********************************************************************
!
! Store this element if it is diagonal or its value is greater than
! CUTOFF
!
         IF ( (IR == IC) .OR. (DABS (ELEMNT) > CUTOFF) ) THEN
            NELC       = NELC + 1
            EMT(NELC)  = ELEMNT
            IROW(NELC) = IR
         ENDIF
!
   85    CONTINUE
! zou
!        print *, ic,SLF_EN(IC)
         IF(LSE) EMT(NELC) = EMT(NELC) + SLF_EN(IC)
! zou
!
!   This column is done; write it to disk
!
         WRITE (imcdf) NELC, ELSTO, (EMT(IR), IR = 1, NELC),       &
                                   (IROW(IR), IR = 1, NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This EAV (and the above EMT) does not have ELSTO.
         EAV = EAV + EMT(NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Output IC on screen to show how far the calculation has proceeded
!
!-----------------------------------------------------------------------
!cjb      if (mod(ic-1,100).eq.0) then
!cjb         PRINT '(A4,I10,A9,I10,A9,I5,A8,I5)', 'Row ', ic,       &
!cjb                 ' nnonz = ', nelc,                             &
!cjb                 ' block = ', jblock,' myid = ', myid
!cjb      ENDIF
!
!cjb          CALL CONVRT (IC,CNUM,LCNUM)
!cjb          if (mod(IC,100).eq.0) then
!cjb             PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
!cjb          end if
!cjb
!
!cjb MPI progress begin --------------------------------------------------
!
!cjb just started
!cychen: report only from myid.eq.0
      if (myid .eq. 0 ) then
       IF (IC .LE. nprocs) then
            PRINT '(A6,I10,A9,I10,A9,I5,A8,I5)', 'Start ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
!
!cjb progress
       IF (IC .GT. nprocs .and. IC .LT. NCF-nprocs .and. &
           MOD (IC-1,100*nprocs) .EQ. myid) then
         if (myid .eq. 0) then
            PRINT '(A6,I10,A9,I10,A9,I5,A8,I5)', 'Done ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         else
            PRINT '(A6,I10,A9,I10,A9,I5,A8,I5)', 'Row  ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         endif
      endif
!
!cjb almost finished
       IF (IC .GT. NCF-nprocs) then

            PRINT '(A7,I10,A9,I10,A9,I5,A8,I5)', 'Finish ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
      endif
!cjb MPI progress end ----------------------------------------------------
!
!
!   Update the counter for the total number of elements
!
         NELMNT = NELMNT + NELC
!
   10 CONTINUE
!***********************************************************************
!
!   Deallocate storage for the arrays in /BUFFER/
!
      CALL ALCBUF (3)

!     ...Locals
      CALL DALLOC (EMT, 'EMT', 'SETHAM')
      CALL DALLOC (IROW, 'IROW', 'SETHAM')

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      CUTOFFtmp = CUTOFF
      NCOEItmp = NCOEI
      NCOECtmp = NCOEC
      NCTEItmp = NCTEI
      NCTECtmp = NCTEC
      NTPItmp = NTPI
      NMCBPtmp = NMCBP
      NCOREtmp = NCORE
      NVPItmp = NVPI
      NKEItmp = NKEI
      NVINTItmp = NVINTI
      NELMNTtmp = NELMNT
      NCFtmp = NCF

      RETURN
      END SUBROUTINE SETHAM
