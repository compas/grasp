!***********************************************************************
      PROGRAM rfinetune
!
!    This program fine-tune the Hamiltonian matrix from rci(rci_mpi) program
!
!    1)Read the Hamiltonian matrix of jj-coupling from rcixxx.res
!    2)Extract the part belonging to the MR
!    3)Transform only this part to LSJ-coupling
!    4)Fine-tune the Hamiltonian matrix of LSJ-coupling
!    5)Transform back to jj-coupling
!    6)Create rcixxx.resnew to replace the fine-tuned part of rcixxx.res
!
!    Variables:
!    H_jj     -- Hamiltonian matrix of jj-coupling
!    H_LSJ    -- Hamiltonian matrix of LSJ-coupling
!    Tmatrix  -- Transformation matrix get from
!                jj2lsj_2022 by Gediminas Gaigalas
!    I_num    -- The number of CSFs in MR
!    ind      -- Number of disk files rcixxx.res, ind = NPROCS
!                                                 Nres= MYID
!    ELSTO2   -- ELSTO in node-0
!                In rci_mpi, only node-0 has the correct, non-zero elsto.
!
!
!    Writen by Yanting Li  30/1/2022
!
!***********************************************************************
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP

      IMPLICIT NONE

      CHARACTER(LEN=6) :: STRING
      CHARACTER(LEN=128) :: NAME
      CHARACTER(LEN=8) :: G92MIX
      CHARACTER :: String2*18
      CHARACTER*300 :: line1,line2
      CHARACTER*500 :: dirstring, DEFNAM, DEFNAMNEW
      CHARACTER :: RECORD*15, S_orbitals*1070
      CHARACTER (LEN = 3) :: idstring
      CHARACTER (LEN = 1) :: transflag

      INTEGER :: I, J, IOS, IR, ind, Nres, nres0
      INTEGER :: IMCDF, IMCDFNEW, NELECR, NCFRES, NWRES, NBLOCKRES, NPARM, N
      INTEGER :: NNUC, NBLOCK, NW
      INTEGER :: NP10, NCF, ICCUT, MYID, NPROCS, NELC, JBLOCK
      INTEGER :: LENNAME, IERR, NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM, JBDUM,&
                 NCFINBLKDUM, NEVINBLKDUM, IATJPDUM, IDUM, LS_number, jj_number
      INTEGER :: ijj, iLS, JB, ii, jj, K, nfine, l, ll, JC, JR, NELC2
      INTEGER :: I_num, Irest
      REAL(DOUBLE) :: wa_transformation, conv, Ediff
      INTEGER, ALLOCATABLE :: ICCUTBLK(:), MF(:), IROW(:)
      INTEGER, ALLOCATABLE :: INELC(:), IROWSAV(:,:)

      DOUBLE PRECISION :: Z, EMN, C, WFACT, RNT, H, HP, ELSTO, ELSTO2
      DOUBLE PRECISION, ALLOCATABLE :: PARM(:), ZZ(:), R(:), RP(:), RPOR(:)
      DOUBLE PRECISION, ALLOCATABLE :: E(:), GAMA(:), PZ(:),EMT(:),EMTNEW(:)
      DOUBLE PRECISION, ALLOCATABLE :: PF(:,:), QF(:,:), H_JJ(:,:), H_LSJ(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Tmatrix(:,:)

      LOGICAL :: LFORDR, LTRANS, LVP, LNMS, LSMS, FOUND

      conv = 2.1947463136320e5
      OPEN(UNIT=215,FILE='rfinetune.log',STATUS='UNKNOWN')

      WRITE (*,*) '****************************************************'
      WRITE (*,*) 'RFINETUNE'
      WRITE (*,*) 'This is the rfinetune program'
      WRITE (*,*) 'Input files:  rci.res, name.lsj.T, name.lsj.c,name.c'
      WRITE (*,*) 'Output files: rci.resnew'
      WRITE (*,*) '****************************************************'
! Ask about transformation matrix from serial calculation or parallel
      WRITE (*,*) 'Transformation matrix is from calculation of:'
      WRITE (*,*) '0--serial'
      WRITE (*,*) '1--parallel'
      READ (*,*) transflag
      WRITE(215,*) transflag
      IF (transflag .NE. '0' .AND. transflag .NE. '1') THEN
         WRITE(*,*) 'Your input must be "0" or "1". Try again!'
         STOP
      ENDIF
! ---------G.G: Read the transformation file .lsj.T --------------
!
      WRITE(*,*) 'Name of MR state: '
      READ(*,*) NAME
      LENNAME = LEN_TRIM(NAME)
      WRITE(215,'(A)') NAME(1:LENNAME)
      OPEN(59,FILE=name(1:LENNAME)//'.lsj.T',FORM='UNFORMATTED'&
             ,STATUS='OLD',IOSTAT=IERR)
      if (ierr /= 0) then
         print *, 'Error when opening ',name, '.lsj.T'
         stop
      end if
      OPEN(95,FILE=name(1:LENNAME)//'.lsj.c',FORM='FORMATTED'&
             ,STATUS='OLD',IOSTAT=IERR)
      if (ierr /= 0) then
         print *, 'Error when opening ',name, '.lsj.c'
         stop
      end if
!    *   *   *
      READ (59, IOSTAT=IERR) G92MIX
      read (59) NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
      !if(NELECDUM /= NELECR .or. NCFTOTDUM /= NCFRES .or. NWDUM /= NW
      !&
      !        .or. NBLOCKDUM /= NBLOCK) then
      !   print*, NELECR, NCFRES, NW, NBLOCK
      !   print*, NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
      !   print*, "Wrong transformation file *.lsj.T"
      !   close(59)
      !   stop
      !end if
      do JB = 1, NBLOCKDUM
         !IATJPDUM = ITJPO(JB)
         read (59) JBDUM, NCFINBLKDUM, NEVINBLKDUM, IATJPDUM
         if (JBDUM /= JB) then
            print*, JB,JBDUM
            print*, "Wrong transformation file *.lsj.T"
            close(59)
            stop
         end if
      end do
!   The head of name.lsj.c
      read(95,*)
      read(95,*)
!   Read the MR CSFs list name.c, obtain the number of CSFs I_num
      OPEN(21,FILE=name(1:LENNAME)//'.c',FORM='FORMATTED'&
             ,STATUS='OLD',IOSTAT=IERR)
      read(21,'(1A15)')RECORD
      IF ( RECORD(1:15)/='Core subshells:') THEN
         WRITE (*, *) 'Not a Configuration Symmetry List File;'
         CLOSE(21)
      ENDIF
      read(21,'(A)')
      read(21,'(A)')
      read(21,'(A)')
      read(21,'(A)')

! rci.res from parallel calculation...
      IF (transflag == '1') THEN
        WRITE (*,*) 'Transformation matrix is from parallel&
                         calculation...'
!        dirstring = '/home/ytli/tmp/'
        WRITE(*,*) 'Name of the temporary directory: (e.g.',"'",&
                 '/home/user/tmp/',"'",')'
        READ(*,*) dirstring
        WRITE(215,'(a1, A, a1)') "'", trim(dirstring), "'"
! ------Inquire and read rcixxx.res----
        ind = 0
        DO
          WRITE (idstring, '(I3.3)') ind
          DEFNAM = trim(dirstring)// idstring //&
                   '/rci' // idstring // '.res'
          DEFNAMNEW = trim(dirstring)// idstring //&
                   '/rci' // idstring // '.resnew'
          INQUIRE(FILE=DEFNAM, EXIST=FOUND)
          IF ( FOUND ) THEN
            IMCDF = 99 + ind
            IMCDFNEW = 599 + ind
!   Open rcixxx.res and rcixxx.resnew
            OPEN(UNIT=IMCDF,FILE=DEFNAM,FORM='UNFORMATTED',STATUS='OLD')
            OPEN(UNIT=IMCDFNEW,FILE=DEFNAMNEW,FORM='UNFORMATTED',&
                  STATUS='UNKNOWN')
! --------- THESE READ STATEMENTS ARE FROM LODRES --------------
!  Read the header and check correctness
            READ(IMCDF) STRING
            WRITE(IMCDFNEW) STRING
            IF (STRING.NE.'R92RES') THEN
               WRITE(*,*) 'Not an rci.res file'
               STOP
            END IF
!   Read the basic parameters of the electron cloud
            READ(IMCDF) NELECR, NCFRES, NWRES, NBLOCKRES
            WRITE(IMCDFNEW) NELECR, NCFRES, NWRES, NBLOCKRES
            NBLOCK = NBLOCKRES   ! If a complete run then these are equal
            NW = NWRES           !
!   Read the nuclear parameters
            READ(IMCDF) Z, EMN
            WRITE(IMCDFNEW) Z, EMN
            ALLOCATE(PARM(10))
            READ(IMCDF) NPARM, (PARM(I),I=1,NPARM)
            WRITE(IMCDFNEW) NPARM, (PARM(I),I=1,NPARM)
            ALLOCATE(ZZ(10000))
            READ(IMCDF) N, (ZZ(I),I=1,N), NNUC
            WRITE(IMCDFNEW) N, (ZZ(I),I=1,N), NNUC
            ALLOCATE(ICCUTBLK(NBLOCK))
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
            READ(IMCDF) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, &
              WFACT, LVP, LNMS, LSMS
            WRITE(IMCDFNEW) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS,&
              WFACT, LVP, LNMS, LSMS
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
            NP10 = N + 10
            ALLOCATE(R(NP10))
            ALLOCATE(RP(NP10))
            ALLOCATE(RPOR(NP10))
            READ(IMCDF) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), &
                    (RPOR(I),I=1,NP10)
            WRITE(IMCDFNEW) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), &
                    (RPOR(I),I=1,NP10)
!   Read the orbital wavefunctions and the associated arrays
            ALLOCATE(E(NW))
            ALLOCATE(GAMA(NW))
            ALLOCATE(PZ(NW))
            ALLOCATE(MF(NW))
            ALLOCATE(PF(NNNP,NW))
            ALLOCATE(QF(NNNP,NW))
            DO J = 1, NW
               READ(IMCDF) E(J), GAMA(J), PZ(J), MF(J)
               WRITE(IMCDFNEW) E(J), GAMA(J), PZ(J), MF(J)
               READ(IMCDF) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
               WRITE(IMCDFNEW) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
            END DO
            DEALLOCATE(PARM)
            DEALLOCATE(ZZ)
            DEALLOCATE(ICCUTBLK)
            DEALLOCATE(R)
            DEALLOCATE(RP)
            DEALLOCATE(RPOR)
            DEALLOCATE(E)
            DEALLOCATE(GAMA)
            DEALLOCATE(PZ)
            DEALLOCATE(MF)
            DEALLOCATE(PF)
            DEALLOCATE(QF)
            ind = ind + 1
          ELSE
            GOTO 73
          ENDIF
        ENDDO
!
!   Print on screen
!
        73 WRITE(*,*) 'There are ', ind, 'files in the temp &
          directory...'
!
!  Loop over block
        DO JBLOCK = 1,NBLOCK
           WRITE(*,*) 'BLOCK',JBLOCK
           I_num = 0
           DO
             read(21,'(A)', end=16) S_orbitals
             IF (S_orbitals(1:2) == ' *') GOTO 16
             I_num = I_num + 1
             read(21,'(A)') S_orbitals
             read(21,'(A)') S_orbitals
           ENDDO

!  Allocate for fine-tuning matrix
        16 ALLOCATE(EMTNEW(I_num*I_num))
           ALLOCATE(H_JJ(I_num,I_num))
           ALLOCATE(H_LSJ(I_num,I_num))
           ALLOCATE(Tmatrix(I_num,I_num))
           ALLOCATE(IROWSAV(I_num,I_num))
           ALLOCATE(INELC(I_num))

           H_JJ = 0.d0
           H_LSJ = 0.d0
           Tmatrix = 0.d0
!
           read(59) String2, IDUM
           if ( String2(1:18) /= ' *   Block Number=' .or. IDUM /= jblock) then
              print*, "Error in transformation file *.lsj.T"
              stop
           end if

           do LS_number = 1, I_num
              do jj_number = 1, I_num
                 read(59) ijj, iLS, wa_transformation
                 if(ijj /= jj_number .or. iLS /= LS_number ) then
                    print*, "Error in jj2lsj transformation file"
                    stop
                 end if
                 Tmatrix(iLS,ijj) = wa_transformation
              end do
           end do
! Loop over rcixxx.res to read
           DO Nres = 0, ind-1
              IMCDF = 99 + Nres
              IMCDFNEW = 599 + Nres
              READ(IMCDF, IOSTAT=IOS) NCF, ICCUT, MYID, NPROCS
              IF (IOS /= 0) THEN
                WRITE(*,*)'READ WRONG IN rci.res, myid=', Nres
              ENDIF
              WRITE(IMCDFNEW) NCF, ICCUT, MYID, NPROCS
           ENDDO
!
           ALLOCATE(EMT(2*NCF))
           ALLOCATE(IROW(2*NCF))
! Loop over rcixxx.res to obtain H_JJ
           INELC(:) = 0
           IROWSAV(:,:) = 0
           DO Nres = 0, ind-1
             IF ( Nres .GE. I_num) GOTO 776
             IMCDF = 99 + Nres
             IMCDFNEW = 599 + Nres
             DO I = Nres + 1, I_num, NPROCS
                READ(IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), &
                   (IROW(IR),IR=1,NELC)
                INELC(I) = NELC
                IF ( Nres == 0 ) ELSTO2 = ELSTO
!  Transfer to dense matrix H_JJ
                DO IR = 1, NELC - 1
                   H_JJ(I,IROW(IR)) = EMT(IR)
                   H_JJ(IROW(IR),I) = EMT(IR)
                   IROWSAV(IR,I) = IROW(IR)
                END DO
                H_JJ(I,IROW(NELC)) = EMT(NELC) + ELSTO2
                IROWSAV(NELC,I) = IROW(NELC)
             END DO
           ENDDO
! ------- CREATE H_LSJ ----------------------
           776 H_LSJ = matmul(matmul(Tmatrix, H_JJ),transpose(Tmatrix))
           !write(*,*)'H_LSJ matrix........'
           !Do ii = 1, I_num
           !  write(388,*)ii, H_LSJ(ii,:)
           !ENDDO
! Read LSJ-coulping CSFs
           DO jj = 1,I_num
             read(95,'(a)') line1
             if (line1(1:2).eq.' *') then
               read(95,'(a)') line1
             endif
             read(95,'(a)') line2
             write(*,*) trim(line1)
             write(*,*) trim(line2)
           ENDDO

! ---------------  FINE-TUNE  ---------------------------
           67 write(*,*) 'How many diagonal elements should &
                         be fine-tuned:'
           read(*,*) nfine
           WRITE(215,*) nfine
           if ( nfine .GT. I_num ) then
              write(*,*) 'Please enter a number less than or equal to'&
                       , I_num
              goto 67
           endif
!
           DO K = 1, nfine
             write(*,*)'Give the serial number of the CSF in LSJ-couping &
               you should fine-tune together with the energy change &
               in cm-1'
             read(*,*)l,Ediff
             WRITE(215,'(I3, a1, f9.2)') l, ",", Ediff
             H_LSJ(l,l) = H_LSJ(l,l) + Ediff/conv
           ENDDO
!
           !write(*,*)'H_LSJ matrix after fintuning........'
           !Do J = 1, I_num
           !  write(389,*)j, H_LSJ(j,:)
           !ENDDO
!
! ------------- TRANFORMN BACK TO H_JJ  ------
           H_JJ = matmul(matmul(transpose(Tmatrix), H_LSJ),Tmatrix)
           !write(*,*)'H_JJ matrix after finetuning........'
           !Do ll = 1, I_num
           !  write(*,*)ll, H_JJ(ll,:)
           !ENDDO

! -------- SUBTRACT ELSTO AND WRITE TO RCI.RESNEW ------
           DO Nres = 0, ind-1
            IMCDF = 99 + Nres
            IMCDFNEW = 599 + Nres
! nres0--count MR for each file
            nres0 = 0
            EMTNEW(:) = 0.d0
            IF ( Nres .LT. I_num) THEN
              DO JC = Nres + 1, I_num, NPROCS
                NELC2 = INELC(JC)
                DO JR = 1,NELC2-1
                   EMTNEW(JR) = H_JJ(JC,IROWSAV(JR,JC))
                ENDDO
                EMTNEW(NELC2) = H_JJ(JC,IROWSAV(NELC2,JC)) - ELSTO2

                WRITE (IMCDFNEW) NELC2, ELSTO2, (EMTNEW(IR), IR = 1, NELC2),&
                                         (IROWSAV(IR,JC), IR = 1, NELC2)
                nres0 = nres0 + 1
              ENDDO
!  Copy the rest matrixs from rci.res to rci.resnew
              DO Irest = nres0 * ind + Nres + 1, NCF, NPROCS
                 READ(IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), &
                    (IROW(IR),IR=1,NELC)
                 WRITE (IMCDFNEW) NELC, ELSTO, (EMT(IR), IR = 1, NELC),&
                                        (IROW(IR), IR = 1, NELC)
              ENDDO
            ELSE
!  if Nres >= I_num, there is no fine-tuning for this node
!  just copy them directly
             DO Irest = Nres + 1, NCF, NPROCS
                READ(IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), &
                   (IROW(IR),IR=1,NELC)
                WRITE (IMCDFNEW) NELC, ELSTO, (EMT(IR), IR = 1, NELC),&
                                       (IROW(IR), IR = 1, NELC)
             END DO
            ENDIF
           ENDDO

!  Deallocate EMT and IROW
           DEALLOCATE(EMT)
           DEALLOCATE(EMTNEW)
           DEALLOCATE(IROW)
           DEALLOCATE(Tmatrix)
           DEALLOCATE(H_LSJ)
           DEALLOCATE(H_JJ)
           DEALLOCATE(IROWSAV)
           DEALLOCATE(INELC)

        END DO

        close(21)
        close(59)
        close(95)
        close(215)
        DO Nres = 0, ind-1
           IMCDF = 99 + Nres
           IMCDFNEW = 599 + Nres
           close(IMCDF)
           close(IMCDFNEW)
        ENDDO

        write(*,*)'Created rcixxx.resnew in ', trim(dirstring)

!  rci.res from serial calculation ...
      ELSEIF (transflag == '0') THEN
        WRITE (*,*) 'Transformation matrix is from serial &
                         calculation...'
        IMCDF = 26
        IMCDFNEW = 36

        OPEN(UNIT=IMCDF,FILE='rci.res',FORM='UNFORMATTED',STATUS='OLD')
        OPEN(UNIT=IMCDFNEW,FILE='rci.resnew',FORM='UNFORMATTED',STATUS=&
                           'UNKNOWN')

! Read the header and check correctness
!
        READ(IMCDF) STRING
        WRITE(IMCDFNEW) STRING
        IF (STRING.NE.'R92RES') THEN
           WRITE(*,*) 'Not an rci.res file'
           STOP
        END IF
!
!   Read the basic parameters of the electron cloud
!
        READ(IMCDF) NELECR, NCFRES, NWRES, NBLOCKRES
        WRITE(IMCDFNEW) NELECR, NCFRES, NWRES, NBLOCKRES
        NBLOCK = NBLOCKRES   ! If a complete run then these are equal
        NW = NWRES           !
! ------------------------------------------------------------
!   Read the nuclear parameters
!
        READ(IMCDF) Z, EMN
        WRITE(IMCDFNEW) Z, EMN
        !WRITE(*,*) Z, EMN

        ALLOCATE(PARM(10))
        READ(IMCDF) NPARM, (PARM(I),I=1,NPARM)
        WRITE(IMCDFNEW) NPARM, (PARM(I),I=1,NPARM)
        !WRITE(*,*) NPARM, (PARM(I),I=1,NPARM)

        ALLOCATE(ZZ(10000))
        READ(IMCDF) N, (ZZ(I),I=1,N), NNUC
        WRITE(IMCDFNEW) N, (ZZ(I),I=1,N), NNUC
        !WRITE(*,*) N, (ZZ(I),I=1,N), NNUC

!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
        ALLOCATE(ICCUTBLK(NBLOCK))
        READ(IMCDF) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, WFACT, LVP, &
           LNMS, LSMS
        WRITE(IMCDFNEW) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, WFACT, LVP, &
           LNMS, LSMS
        !WRITE(*,*) C, LFORDR, (ICCUTBLK(I),I=1,NBLOCK), LTRANS, WFACT,
        !LVP, &
        !   LNMS, LSMS
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
        NP10 = N + 10
        ALLOCATE(R(NP10))
        ALLOCATE(RP(NP10))
        ALLOCATE(RPOR(NP10))
        READ(IMCDF) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), &
                (RPOR(I),I=1,NP10)
        WRITE(IMCDFNEW) RNT, H, HP, (R(I),I=1,NP10), (RP(I),I=1,NP10), &
                (RPOR(I),I=1,NP10)

!   Read the orbital wavefunctions and the associated arrays
!

        ALLOCATE(E(NW))
        ALLOCATE(GAMA(NW))
        ALLOCATE(PZ(NW))
        ALLOCATE(MF(NW))
        ALLOCATE(PF(NNNP,NW))
        ALLOCATE(QF(NNNP,NW))
        DO J = 1, NW
           READ(IMCDF) E(J), GAMA(J), PZ(J), MF(J)
           WRITE(IMCDFNEW) E(J), GAMA(J), PZ(J), MF(J)
           !WRITE(*,*) E(J), GAMA(J), PZ(J), MF(J)
           READ(IMCDF) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
           WRITE(IMCDFNEW) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
           !WRITE(*,*) (PF(I,J),I=1,MF(J)), (QF(I,J),I=1,MF(J))
        END DO
!
        DO JBLOCK = 1,NBLOCK
           WRITE(*,*) 'BLOCK',JBLOCK
! --------- THESE READ STATEMENTS ARE FROM GENMAT --------------
           READ(IMCDF) NCF, ICCUT, MYID, NPROCS
           WRITE(IMCDFNEW) NCF, ICCUT, MYID, NPROCS
           I_num = 0
           DO
             read(21,'(A)', end=17) S_orbitals
             IF (S_orbitals(1:2) == ' *') GOTO 17
             I_num = I_num + 1
             read(21,'(A)') S_orbitals
             read(21,'(A)') S_orbitals
           ENDDO
           write(*,*)'I_num=', I_num

!   We allocate NCF x NCF for EMT and IROW which is more than enough
           17 ALLOCATE(EMT(NCF*NCF))
           ALLOCATE(EMTNEW(I_num*I_num))
           ALLOCATE(IROW(NCF*NCF))
           ALLOCATE(H_JJ(I_num,I_num))
           ALLOCATE(H_LSJ(I_num,I_num))
           ALLOCATE(Tmatrix(I_num,I_num))
           ALLOCATE(IROWSAV(I_num*I_num,I_num))
           ALLOCATE(INELC(I_num))

           H_JJ = 0.d0
           H_LSJ = 0.d0
           Tmatrix = 0.d0
!
           read(59) String2, IDUM
           if ( String2(1:18) /= ' *   Block Number=' .or. IDUM /= jblock) then
              print*, "Error in transformation file *.lsj.T"
           end if


           do LS_number = 1, I_num
              do jj_number = 1, I_num
                 read(59) ijj, iLS, wa_transformation
                 if(ijj /= jj_number .or. iLS /= LS_number ) then
                    print*, "Error in jj2lsj transformation file"
                    stop
                 end if
                 Tmatrix(iLS,ijj) = wa_transformation
              end do
           end do
!
!
           INELC(:) = 0
           IROWSAV(:,:) = 0
           DO I = MYID + 1, I_num, NPROCS
              READ(IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), &
                 (IROW(IR),IR=1,NELC)
              INELC(I) = NELC
              DO k = 1, NELC
                 IROWSAV(K,I) = IROW(k)
              ENDDO

!  Transfer to dense matrix H_JJ

              DO IR = 1, NELC - 1
                 H_JJ(I,IROW(IR)) = EMT(IR)
                 H_JJ(IROW(IR),I) = EMT(IR)
              END DO
              H_JJ(I,IROW(NELC)) = EMT(NELC) + ELSTO  !
           END DO

! ------- CREATE H_LSJ ----------------------
           H_LSJ = matmul(matmul(Tmatrix, H_JJ),transpose(Tmatrix))
           !write(*,*)'H_LSJ matrix........'
           !Do ii = 1, I_num
           !  write(*,*)ii, H_LSJ(ii,:)
           !ENDDO
! Read LSJ csl file
           DO jj = 1,I_num
             read(95,'(a)') line1
             if (line1(1:2).eq.' *') then
               read(95,'(a)') line1
             endif
             read(95,'(a)') line2
             write(*,*) 'No. in LSJ-couping =', jj
             print*, trim(line1)
             print*, trim(line2)
           ENDDO

! ---------------  FINE TUNE  ---------------------------
           68 write(*,*) 'How many diagnoal elements should be &
                       fine-tuned:'
           read(*,*) nfine
           WRITE(215,*) nfine
           if ( nfine .GT. I_num ) then
              write(*,*) 'Please enter a number less than or equal to'&
                       , I_num
              goto 68
           endif
!
           DO K = 1, nfine
            write(*,*)'Give the serial number of the CSF in LSJ-couping you &
             should fine-tune together with the energy change in cm-1'
            read(*,*)l,Ediff
            WRITE(215,'(I3 a1, f9.2)') l, ",", Ediff
            H_LSJ(l,l) = H_LSJ(l,l) + Ediff/conv
           ENDDO
!
           !write(*,*)'H_LSJ matrix after fintuning........'
           !Do J = 1, I_num
           !  write(*,*)j, H_LSJ(j,:)
           !ENDDO
!
! ------------- TRANFORMN BACK TO H_JJ  ------
           H_JJ = matmul(matmul(transpose(Tmatrix), H_LSJ),Tmatrix)
           !write(*,*)'H_JJ matrix after finetuning........'
           !Do ll = 1, I_num
           !  write(*,*)ll, H_JJ(ll,:)
           !ENDDO

! -------- SUBTRACT ELSTO AND WRITE TO RCI.RESNEW ------
           EMTNEW(:) = 0.d0
           DO JC = MYID + 1, I_num, NPROCS
              NELC2  = INELC(JC)
              DO JR = 1,NELC2-1
                 EMTNEW(JR) = H_JJ(JC,IROWSAV(JR,JC))
              ENDDO
              EMTNEW(NELC2) = H_JJ(JC,IROWSAV(NELC2,JC)) - ELSTO
              WRITE (imcdfnew) NELC2, ELSTO, (EMTNEW(IR), IR = 1, NELC2),&
                                       (IROWSAV(IR,JC), IR = 1, NELC2)
           ENDDO
!  Copy the rest of matrixs from rci.res to rci.resnew
           DO Irest = MYID + I_num + 1, NCF, NPROCS
              READ(IMCDF) NELC, ELSTO, (EMT(IR),IR=1,NELC), &
                 (IROW(IR),IR=1,NELC)
              WRITE (imcdfnew) NELC, ELSTO, (EMT(IR), IR = 1, NELC),&
                                     (IROW(IR), IR = 1, NELC)
           END DO


!  Deallocate EMT and IROW
           DEALLOCATE(EMT)
           DEALLOCATE(EMTNEW)
           DEALLOCATE(IROW)
           DEALLOCATE(Tmatrix)
           DEALLOCATE(H_LSJ)
           DEALLOCATE(H_JJ)
           DEALLOCATE(IROWSAV)
           DEALLOCATE(INELC)

        END DO
        close(IMCDF)
        close(IMCDFNEW)
        close(21)
        close(59)
        close(95)
        close(215)

        write(*,*)'Created rci.resnew'
      ENDIF

      END PROGRAM
