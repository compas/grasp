!***********************************************************************
!                                                                      *
      SUBROUTINE QQSORT(NFILE, NUMBER, KSTART, NAME, KAMAX) 
!                                                                      *
!     The list of unique integrals (j,i) is formed in the order of     *
!     increasing symmetry, i.e. with j .le. i.                         *
!     With each integral (INTGRL) there is a pointer (INTPTR)          *
!     that points to the last element in the array of coefficients     *
!     for that integral.                                               *
!                                                                      *
!     On Exit                                                          *
!                                                                      *
!     INTGRL -- array of integrals (packed form of j,i)                *
!     NINTG  -- number of integrals                                    *
!     INTPTR -- array of pointers to last element in list of coeff.    *
!     CNN    -- array of coefficients                                  *
!     NCOEFF -- number of coefficients                                 *
!     JANN   -- row of matrix                                          *
!     JBNN   -- column of matrix                                       *
!                                                                      *
!     Written by Per Jonsson                                           *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:38:04   1/ 6/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: KEYORB
      USE memory_man
      USE orb_C
      USE default_C
      USE mcpdata_C
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NFILE 
      INTEGER, INTENT(IN) :: NUMBER 
      INTEGER  :: KSTART 
      INTEGER, INTENT(IN) :: KAMAX 
      CHARACTER  :: NAME*24 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
      INTEGER, PARAMETER :: NF = 200 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, L, IR, JA, JB, INT, J, INTT
      REAL(DOUBLE) :: CN 
!-----------------------------------------------
!
      REWIND (NFILE + 80) 
 
      NCOEFF = NUMBER 
!
!   Sort the list
!
!   Allocate storage for all required arrays these arrays are then deallocated
!   in mcp
!
      CALL ALLOC (JANN,   NCOEFF, 'JANN',   'QQSORT') 
      CALL ALLOC (JBNN,   NCOEFF, 'JBNN',   'QQSORT') 
      CALL ALLOC (INTGRL, NCOEFF, 'INTGRL', 'QQSORT') 
      CALL ALLOC (CNN,    NCOEFF, 'CNN',    'QQSORT') 
      CALL ALLOC (INTPTR, NCOEFF, 'INTPTR', 'QQSORT') 
!
!   Read arrays into memory from NFILE
!
      DO I = 1, NCOEFF 
         READ (NFILE + 80) JANN(I), JBNN(I), INTGRL(I), CNN(I) 
      END DO 
!
!   Sort INTGRL into ascending order using the heapsort algorithm;
!   (Numerical recepies page 231.) move the associated members of
!   of JANN and JBNN in the same
!   manner;
!
      IF (NCOEFF > 1) THEN 
!
         L = NCOEFF/2 + 1 
         IR = NCOEFF 
    2    CONTINUE 
         IF (L > 1) THEN 
            L = L - 1 
            JA = JANN(L) 
            JB = JBNN(L) 
            INT = INTGRL(L) 
            CN = CNN(L) 
         ELSE 
            JA = JANN(IR) 
            JB = JBNN(IR) 
            INT = INTGRL(IR) 
            CN = CNN(IR) 
            JANN(IR) = JANN(1) 
            JBNN(IR) = JBNN(1) 
            INTGRL(IR) = INTGRL(1) 
            CNN(IR) = CNN(1) 
            IR = IR - 1 
            IF (IR == 1) THEN 
               JANN(1) = JA 
               JBNN(1) = JB 
               INTGRL(1) = INT 
               CNN(1) = CN 
               GO TO 4 
            ENDIF 
         ENDIF 
         I = L 
         J = L + L 
    3    CONTINUE 
         IF (J <= IR) THEN 
            IF (J < IR) THEN 
               IF (INTGRL(J) < INTGRL(J+1)) J = J + 1 
            ENDIF 
            IF (INT < INTGRL(J)) THEN 
               JANN(I) = JANN(J) 
               JBNN(I) = JBNN(J) 
               INTGRL(I) = INTGRL(J) 
               CNN(I) = CNN(J) 
               I = J 
               J = J + J 
            ELSE 
               J = IR + 1 
            ENDIF 
            GO TO 3 
         ENDIF 
         JANN(I) = JA 
         JBNN(I) = JB 
         INTGRL(I) = INT 
         CNN(I) = CN 
         GO TO 2 
      ENDIF 
!
 
!   Sorting complete; close the file
!
!ww    4 CLOSE (80+NFILE)
    4 CONTINUE 
      NINTG = 1 
      INTT = INTGRL(1) 
!
      DO I = 1, NCOEFF 
         IF (INTGRL(I) == INTT) CYCLE  
         INTPTR(NINTG) = I - 1 
         NINTG = NINTG + 1 
         INTT = INTGRL(I) 
      END DO 
 
      INTPTR(NINTG) = NCOEFF 
!
      DO I = 1, NINTG 
         INTGRL(I) = INTGRL(INTPTR(I)) 
      END DO 
!
!  If output option is set dump the data on file
!
      IF (NDUMP == 1) THEN 
!
!  If first set of data open the file and print
!  some data to later be able to identify the file
!
!        IF (NFILE.EQ.1.AND.IBLK.EQ.1) THEN
!         print *, ' open sorted ang. file .TB(NF)', NF
!          J = INDEX(NAME,' ')
!          OPEN (UNIT = NF,FILE=NAME(1:J-1)//'.TB',
!     :          STATUS='UNKNOWN',FORM='UNFORMATTED')
!          REWIND (NF)
!        ENDIF
         IF (NFILE == 1) WRITE (NF) NCF,NW,KAMAX,nprocs,myid
!
!  Print out angular data for this kappa
!
         WRITE (NF) NINTG, NCOEFF 
         DO I = 1, NINTG 
            WRITE (NF) INTGRL(I), INTPTR(I) 
         END DO 
         DO I = 1, NCOEFF 
            WRITE (NF) CNN(I), JANN(I), JBNN(I) 
         END DO 
      ENDIF 
      CALL DALLOC (JANN,   'JANN',   'QQSORT') 
      CALL DALLOC (JBNN,   'JBNN',   'QQSORT') 
      CALL DALLOC (INTGRL, 'INTGRL', 'QQSORT') 
      CALL DALLOC (CNN,    'CNN',    'QQSORT') 
      CALL DALLOC (INTPTR, 'INTPTR', 'QQSORT') 
!
!  Has all data been processed? If so close
!  the file
!
!     IF (NFILE.EQ.KAMAX) CLOSE (NF)
 
!  Debug output
!
!      IF (KSTART.EQ.1) THEN
!        NNNN = 120
!      ELSE
!        NNNN = 140
!      ENDIF
!      WRITE(NNNN+NFILE,*) ' Data for '
!      WRITE(NNNN+NFILE,*) NINTG, ' Integrals'
!      DO I = 1,NINTG
!        IA = INTGRL(I)/KEY
!        IB = MOD(INTGRL(I),KEY)
!        WRITE(NNNN+NFILE,302) NP(IA),NH(IA),NP(IB),NH(IB),INTPTR(I)
!      ENDDO
!      WRITE(NNNN+NFILE,*) NCOEFF,' Coefficients'
!      WRITE(NNNN+NFILE,'(F12.8,2I6)')
!     :     (CNN(I),JANN(I),JBNN(I),I=1,NCOEFF)
!
      RETURN  
!
  301 FORMAT(' T_[',1I2,',',1I2,']',' (',1I2,1A2,',',1I2,1A2,') = ',1P,D19.12) 
  302 FORMAT(' (',1I2,1A2,',',1I2,1A2,')',I6) 
      RETURN  
!
      END SUBROUTINE QQSORT 
