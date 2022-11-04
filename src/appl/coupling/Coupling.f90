!***********************************************************************
!
      PROGRAM Coupling
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use Coupling_constants
      use Coupling_structures
      use Coupling_main
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=50) :: log_file
      character(len=24) :: NAME
      character(len=2)  :: note
      character(len=13) :: couplings
      real(kind=dp)     :: thresh
      integer           :: print_level,K
!-----------------------------------------------
      print *, " "
      print *,                                                         &
      "Coupling: Transformation of ASFs from a LS-coupled CSF basis"
      print *,                                                         &
      "          into differete coupled CSF bases      (Fortran 95)"
      print *,                                                         &
      "          (C) (2022)                G. Gaigalas, A. Kramida."
      print *, " "
      print *,                                                         &
      "Input  files: *.lsj.c, *.lsj.j (ATSP (CPC) or GRASP2K types)"
      print *,                                                         &
      "Output files: *.coup*.*.lbl, *.coup*.sum"
      print *, " "
      print *, " "
    1 WRITE (*,*) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE (*,*) 'Names may not start with a blank'
         GOTO 1
      ENDIF
      call get_parameters                                              &
                (couplings,NAME,print_level,note,thresh)
      call open_files(couplings,NAME,print_level)
      call GLCONS
      call FACTRL(32)
      call main(1, couplings, print_level,thresh)
!GG      call store_R_values(iwrite_log, note,couplings)
      call store_P_values(iwrite_log, note,couplings)
      call iseiti
      print *, " "
      write(*,*)'Coupling: Execution complete.'
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine get_parameters(couplings,NAME,print_level,note,thresh)
!                                                                      *
!     This subroutine asks the user for some parameters:               *
!     couplings,                                                       *
!     printing level,                                                  *
!     note (for the R_values, and P_values output),                    *
!     input data file names.                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(in) :: NAME
      character(len=13), intent(out) :: couplings
      character(len=2), intent(out) :: note
      character(len=4)  :: jval
      integer, intent(out)     :: print_level
      real(kind=dp), intent(out) :: thresh
!-----------------------------------------------
      print_level=1
      note='No'
!      write(*,*)'Specify the data source 1-mchf '
!      write(*,*)'mchf'
!      write(*,*)'Specify the print level : 0, 1 or 2'
!      read(*,'(i1)')print_level
      write(*,*)'Default settings ? (Y/N)'
      read(*,'(a1)') note
      if(note(1:1) == 'Y' .or. note(1:1) == 'y') then
         print_level = 0
      else
         print_level = 2
      end if
      note='  '
      write (*,*)                                                      &
       'Specify the number of coupled shells for evaluation (1,2 or 3):'
      read(*,'(i1)') I_Numb
!      write (*,*) '                                       ',I_Numb
      write (*,'(I2)') I_Numb
      symbol_I_Numb = jval(2*I_Numb)
      log_file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.sum'
      open(iwrite_log, file=log_file, status='unknown')
      write(iwrite_log,*) "      -------------------------------"
      write(iwrite_log,*) "      Summry file of program COUPLING"
      write(iwrite_log,*) "      -------------------------------"
      write(iwrite_log,*) ""
      write(iwrite_log,*)                                              &
             "The CSF input file:                ",NAME(1:K-1)//'.lsj.c'
      write(iwrite_log,*)                                              &
             "The mixing coefficient input file: ",NAME(1:K-1)//'.lsj.j'
      write(iwrite_log,*) "There are",I_Numb,                          &
                                 "(in LS-coupling) coupled shells case."
      if (I_Numb == 1) then
         couplings='ynnnnnnnnnynn'
      else if (I_Numb == 2) then
         couplings='yyyynnnnnynyn'
      else if (I_Numb ==3) then
         couplings='ynnnyyyyynnny'
      else
         write(*,*)                                                   &
           'Not proper number of coupled shells for evaluation.'
         stop
      end if
      write(iwrite_log,*) "Specification of couplings:  ",couplings
      write(iwrite_log,*)                                              &
      "(LS JJ LK JK LS3 LSJ3 LK3 JK3 cLSJ3 LScjj jj1 jj2 jj3)"
      write(iwrite_log,                                                &
      '(3X,A1,4(2X,A1),3X,A1,4X,A1,3X,A1,4X,A1,4(4X,A1))')             &
      couplings(1:1),couplings(2:2),couplings(3:3),couplings(4:4),     &
      couplings(5:5),couplings(6:6),couplings(7:7),couplings(8:8),     &
      couplings(9:9),couplings(10:10),couplings(11:11),                &
      couplings(12:12),couplings(13:13)
      write (*,*)                                                      &
              'What is the value below which an eigenvector composition'
      write (*,*) 'is to be neglected for printing?'
      READ *, thresh
      write (*,*) '                                ',thresh
      write(iwrite_log,*) ""
      write(iwrite_log,*)                                              &
                   'The value below which an eigenvector composition'
      write(iwrite_log,*) 'is to be neglected for printing:',thresh
      write(iwrite_log,*) ""
      end subroutine get_parameters
!
!***********************************************************************
!                                                                      *
      subroutine open_files(couplings,NAME,print_level)
!                                                                      *
!     This subroutine opens the files of input and output              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=13), intent(in) :: couplings
      character(len=24), intent(in) :: NAME
      integer, intent(in) :: print_level
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=50) :: suggestions_file
!print_level 1
      character(len=50) :: classifications_file
      character(len=50) :: classifications_evaluation_file
      character(len=50) :: expansions_file
!print_level 2
      character(len=50) :: log_file2
      character(len=50) :: cfg_lsjj_file,cfg_lslk_file,cfg_lsjk_file
      character(len=51) :: cfg_lsls3_file,cfg_lslsj3_file,cfg_lslk3_file
      character(len=51) :: cfg_lsjk3_file,cfg_lsclsj3_file
      character(len=51) :: cfg_lslscjj_file,cfg_lsjj123_file
      integer           :: K
!-----------------------------------------------
!      write(*,*) 'open_files'
!
      K=INDEX(NAME,' ')
      log_file2 = 'log2.coup'
!
      cfg_lsjj_file     ='ASF_expansions_lsJJ.coup'
      cfg_lslk_file     ='ASF_expansions_lslk.coup'
      cfg_lsjk_file     ='ASF_expansions_jjjk.coup'
      cfg_lsls3_file    ='ASF_expansions_lsls3.coup'
      cfg_lslsj3_file   ='ASF_expansions_lslsj3.coup'
      cfg_lslk3_file    ='ASF_expansions_lslk3.coup'
      cfg_lsjk3_file    ='ASF_expansions_lsjk3.coup'
      cfg_lsclsj3_file  ='ASF_expansions_lsclsj3.coup'
      cfg_lslscjj_file  ='ASF_expansions_lslscjj.coup'
      cfg_lsjj123_file  ='ASF_expansions_lsjj.coup'
!
      expansions_file='ASF_expansions.coup'
!
      classifications_file='All_classifications.coup'
      classifications_evaluation_file='States_classif_evaluations.coup'
!
      open(iread_csf, file=NAME(1:K-1)//'.lsj.c',status='old')
      open(iread_coefs, file=NAME(1:K-1)//'.lsj.j',status='old')
      if(print_level > 3) then
         suggestions_file=                                             &
                   NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.Suggest'
         open(iwrite_suggestions,file=suggestions_file,status='unknown')
      end if
      if(couplings(1:1).eq.'y')                                        &
      open(iwrite_DataBase_LS, file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LS.lbl',status='unknown')
      if(couplings(2:2).eq.'y')                                        &
      open(iwrite_DataBase_JJ,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.JJ.lbl',status='unknown')
      if(couplings(3:3).eq.'y')                                        &
      open(iwrite_DataBase_LK,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LK.lbl',status='unknown')
      if(couplings(4:4).eq.'y')                                        &
      open(iwrite_DataBase_JK,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.JK.lbl',status='unknown')
      if(couplings(5:5).eq.'y')                                        &
      open(iwrite_DataBase_LS3,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LS3.lbl',status='unknown')
      if(couplings(6:6).eq.'y')                                        &
      open(iwrite_DataBase_LSJ3,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LSJ3.lbl',status='unknown')
      if(couplings(7:7).eq.'y')                                        &
      open(iwrite_DataBase_LK3,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LK3.lbl',status='unknown')
      if(couplings(8:8).eq.'y')                                        &
      open(iwrite_DataBase_JK3,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.JK3.lbl',status='unknown')
      if(couplings(9:9).eq.'y')                                        &
      open(iwrite_DataBase_cLSJ3,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.cLSJ3.lbl',status='unknown')
      if(couplings(10:10).eq.'y')                                      &
      open(iwrite_DataBase_LScjj,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.LScjj.lbl',status='unknown')
      if(couplings(11:11).eq.'y')                                      &
      open(iwrite_DataBase_jj123,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.jj.lbl',status='unknown')
      if(couplings(12:12).eq.'y')                                      &
      open(iwrite_DataBase_jj123,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.jj.lbl',status='unknown')
      if(couplings(13:13).eq.'y')                                      &
      open(iwrite_DataBase_jj123,file=NAME(1:K-1)//'.coup'//symbol_I_Numb(1:1)//'.jj.lbl',status='unknown')
      if(print_level.gt.0) then
         open(iwrite_expansions, file=expansions_file, status='unknown')
         open(iwrite_classifications, file=classifications_file,       &
                                                       status='unknown')
         open(iwrite_classifications_evaluations,                      &
                 file=classifications_evaluation_file, status='unknown')
         if(print_level.gt.1) then
            if(couplings(2:2).eq.'y')                                  &
            open(iwrite_cfg_expansions_JJ, file=cfg_lsjj_file,         &
                                                       status='unknown')
            if(couplings(3:3).eq.'y')                                  &
            open(iwrite_cfg_expansions_LK, file=cfg_lslk_file,         &
                                                       status='unknown')
            if(couplings(4:4).eq.'y')                                  &
            open(iwrite_cfg_expansions_JK, file=cfg_lsjk_file,         &
                                                       status='unknown')
            if(couplings(5:5).eq.'y')                                  &
            open(iwrite_cfg_expansions_LS3, file=cfg_lsls3_file,       &
                                                       status='unknown')
            if(couplings(6:6).eq.'y')                                  &
            open(iwrite_cfg_expansions_LSJ3, file=cfg_lslsj3_file,     &
                                                       status='unknown')
            if(couplings(7:7).eq.'y')                                  &
            open(iwrite_cfg_expansions_LK3, file=cfg_lslk3_file,       &
                                                       status='unknown')
            if(couplings(8:8).eq.'y')                                  &
            open(iwrite_cfg_expansions_JK3, file=cfg_lsjk3_file,       &
                                                       status='unknown')
            if(couplings(9:9).eq.'y')                                  &
            open(iwrite_cfg_expansions_cLSJ3, file=cfg_lsclsj3_file,   &
                                                       status='unknown')
            if(couplings(10:10).eq.'y')                                &
            open(iwrite_cfg_expansions_LScjj, file=cfg_lslscjj_file,   &
                                                       status='unknown')
            if(couplings(11:11).eq.'y')                                &
            open(iwrite_cfg_expansions_jj123, file=cfg_lsjj123_file,   &
                                                       status='unknown')
            if(couplings(12:12).eq.'y')                                &
            open(iwrite_cfg_expansions_jj123, file=cfg_lsjj123_file,   &
                                                       status='unknown')
            if(couplings(13:13).eq.'y')                                &
            open(iwrite_cfg_expansions_jj123, file=cfg_lsjj123_file,   &
                                                       status='unknown')
         end if
      end if
!      write(*,*) 'end open_files'
      end subroutine open_files
!
      end program Coupling
