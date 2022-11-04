!
!***********************************************************************
!                                                                      *
      module Coupling_getdata_mchf
!                                                                      *
!     Written by G. Gaigalas,                                          *
!                                                                      *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use Coupling_constants
      use Coupling_structures
      use Coupling_data
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: get_mchf_data
      private :: analy1GG
      private :: analy1
      private :: cfgo2
      private :: cfgo3
      private :: analy_j
      private :: readcoeffs
      private :: print_csf
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      type(set_of_csfs_LS) :: csfs_LS
      integer, dimension(:),pointer:: Iselect
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
      subroutine get_mchf_data(icase, icoupling_nr_in_expansions,     &
                               print_level, nr_of_j)
!                                                                      *
!     This managerial subroutine                                       *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: icase, icoupling_nr_in_expansions, print_level
      integer, intent(inout) :: nr_of_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I_CSF, J
      character(len=72), save :: cfg_line
!-----------------------------------------------
!      write(*,*)'   subroutine get_mchf_data'
!      write(*,*)'Specify csf.inp format:'
!      write(*,*)'                      ATSP2K CPC or GRASP jj2lsj (0/1)'
!      read(*,'(i1)') I_CSF
      I_CSF = 1
!      open(1000+iread_csf, file='Coup.tmp',status='unknown')
!      open(1000+iread_coefs, file='Coup2.tmp',status='unknown')
      open(1000+iread_csf, status='scratch')
      open(1000+iread_coefs, status='scratch')
      if(icase == 1) call analy1GG(nr_of_j,I_CSF)
      call analy1(icase,I_CSF)
      if(csfs_LS%nr_of_csfs.gt.0) then
         allocate(csfs_LS%csfs(csfs_LS%nr_of_csfs))
         if(icase == 1) then
            write(*,*) 'Specify shells for recoupling (no more than 12)'
            read(*,'(a56)') cfg_line
            write(*,*)                                 &
                      "List of the shells included in the recoupling:"
            write(*,*) cfg_line
            write(*,*) ""
            write(iwrite_log,*)                                 &
                      "List of the shells included in the recoupling:"
            write(iwrite_log,*) cfg_line
         end if
         write(*,*) ""
         write(*,*) ""
         write(*,'(A15,I5)') "Symmetry block:",icase
         write(iwrite_log,*)""
         write(iwrite_log,*)""
         write(iwrite_log,'(1x,a15,i5)')"Symmetry block:",icase
         write(iwrite_log,'(1x,a39,i9,a30)')                         &
                   'Total number of CSF in the input file =',        &
                   csfs_LS%nr_of_csfs, ' (number of csf in expansion)'
         write(iwrite_log,*)                                 &
       "---------------------------------------------------------------"
         write(iwrite_log,*)"List of configurations which",  &
            " are eliminated from the calculation"
         write(iwrite_log,*)''
         write(iwrite_log,*)"   No.           Configurations"
         write(iwrite_log,*)                                 &
       "---------------------------------------------------------------"
         if (I_Numb == 1) then
            call cfgo1(I_CSF,cfg_line)
         else if (I_Numb == 2) then
            call cfgo2(I_CSF,cfg_line)
         else if (I_Numb == 3) then
            call cfgo3(I_CSF,cfg_line)
         else
           write(*,*) "To many shell for evaluation"
           stop
         end if
         call analy_j(icase,I_CSF,icoupling_nr_in_expansions)
         if(states%nr_of_states.gt.0 .and.                  &
         all_expansions%coupling_expansions                 &
         (icoupling_nr_in_expansions)%nr_of_expansions.gt.0) then
            allocate(states%states(states%nr_of_states))
            allocate(all_expansions%                                  &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     expansions(all_expansions%                       &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     nr_of_expansions))
         end if
         call readcoeffs(icase,I_CSF,nr_of_j,icoupling_nr_in_expansions)
      end if
      if(print_level.gt.1) call print_csf
      call mchf_iseiti
      close(1000+iread_csf)
      close(1000+iread_coefs)
!      write(*,*)'   end subroutine get_mchf_data'
      contains
!
!***********************************************************************
!                                                                      *
         subroutine mchf_iseiti
!                                                                      *
!***********************************************************************
         implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
         integer :: error
         integer :: icsf
!-----------------------------------------------
         if(associated(Iselect)) deallocate(Iselect, STAT = error)
         if(associated(csfs_LS%csfs)) then
            do icsf=1, csfs_LS%nr_of_csfs
               if(associated(csfs_LS%csfs(icsf)%subc))              &
               deallocate(csfs_LS%csfs(icsf)%subc, STAT = error)
               if(associated(csfs_LS%csfs(icsf)%subc_cfg))          &
               deallocate(csfs_LS%csfs(icsf)%subc_cfg, STAT = error)
               if(associated(csfs_LS%csfs(icsf)%iM1))               &
               deallocate(csfs_LS%csfs(icsf)%iM1, STAT = error)
               if(associated(csfs_LS%csfs(icsf)%iM2))               &
               deallocate(csfs_LS%csfs(icsf)%iM2, STAT = error)
               if(associated(csfs_LS%csfs(icsf)%iJ))               &
               deallocate(csfs_LS%csfs(icsf)%iJ, STAT = error)
            end do
            deallocate(csfs_LS%csfs, STAT = error)
         end if
         end subroutine mchf_iseiti
!
      end subroutine get_mchf_data
!
!***********************************************************************
!                                                                      *
      subroutine analy1GG(nr_of_j,I_CSF)
!                                                                      *
!     This subroutine defines:                                         *
!     csfs_LS%nr_of_csfs  : number of csf's in mchf expansion          *
!     (to be used in definition of array of csf's                      *
!     in subroutine cfgo1)                                             *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: I_CSF
      integer, intent(out) :: nr_of_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=72) :: line
      character(len=55) :: cfg_line
      integer      :: ilinecount, eof
!-----------------------------------------------
!      write(*,*) '   subroutine analy1GG'
      read(iread_csf,'(a72)') line
      read(iread_csf,'(a72)') line
      ilinecount=0
      nr_of_j = 0
      do
         read (iread_csf,'(a55)',iostat=eof) cfg_line
         if (eof /= 0) exit
!GG         if(cfg_line(2:2).eq.'*')read (iread_csf,'(a55)',iostat=eof) cfg_line
         if(cfg_line(1:1).eq.'*' .or. cfg_line(2:2).eq.'*') then
            nr_of_j = nr_of_j +1
            write (2000+iread_csf,'(I55)') ilinecount
            if (eof /= 0) exit
            ilinecount=0
            cycle
         end if
         if(cfg_line(3:3).eq.' '.or.cfg_line(4:4).eq.' ') exit
!         if(cfg_line(3:3).eq.' '.or.cfg_line(4:4).eq.' ') THEN
!           print*, "STOP analy1GG: A"
!           print*, cfg_line
!           STOP
!         END IF
!         if(cfg_line(1:1).ne.' '.or.cfg_line(9:9).ne.' ') exit
         if(cfg_line(1:1).ne.' '.or.cfg_line(9:9).ne.' ') THEN
           print*, "STOP analy1GG: B"
           STOP
         END IF
         read (iread_csf,'(a72)') line
         write (1000+iread_csf,'(a55)') cfg_line
         write (1000+iread_csf,'(a72)') line
         if(I_CSF == 0) THEN
            if(line(1:5).ne.'    '.or.line(6:6).eq.' ') exit
            if(line(7:7).eq.' '.or.line(8:8).eq.' ') exit
         else if(I_CSF == 1) THEN
!            if(line(1:1).ne.' '.or.line(2:2).eq.' ') exit
            if(line(1:1).ne.' '.or.line(2:2).eq.' ') THEN
           print*, "STOP analy1GG: E"
           STOP
         END IF
!            if(line(3:3).eq.' '.or.line(4:4).eq.' ') exit
            if(line(3:3).eq.' '.or.line(4:4).eq.' ') THEN
           print*, "STOP analy1GG: F"
           STOP
         END IF
         else
            write(*,*) ' error in cfg.inp format I_CSF',I_CSF
            stop
         end if
         ilinecount = ilinecount + 1
         if (ilinecount == 0) then
            write(*,*) ' error in cfg.inp format'
            stop
         end if
      end do
      rewind(iread_csf)
      rewind(2000+iread_csf)
!      csfs_LS%nr_of_csfs=ilinecount
!      print*,"analy1 csfs_LS%nr_of_csfs=",csfs_LS%nr_of_csfs
!      write(iwrite_log,'(1x,a39,i9,a30)')                            &
!        'Total number of CSF in the input file =',                   &
!         csfs_LS%nr_of_csfs, ' (number of csf in expansion)'
!      write(*,*) '   end subroutine analy1'
      end subroutine analy1GG
!
!***********************************************************************
!                                                                      *
      subroutine analy1(icase,I_CSF)
!                                                                      *
!     This subroutine defines:                                         *
!     csfs_LS%nr_of_csfs  : number of csf's in mchf expansion          *
!     (to be used in definition of array of csf's                      *
!     in subroutine cfgo1)                                             *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: icase, I_CSF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=72) :: line
      character(len=55) :: cfg_line
      integer      :: ilinecount, eof
!-----------------------------------------------
!      write(*,*) '   subroutine analy1'
      read (2000+iread_csf,'(I55)') csfs_LS%nr_of_csfs
      if(icase == 1) read(iread_csf,'(a72)') line
      if(icase == 1) read(iread_csf,'(a72)') line
      ilinecount=0
      do
         read (iread_csf,'(a55)',iostat=eof) cfg_line
!GG         if(cfg_line(2:2).eq.'*')read (iread_csf,'(a55)',iostat=eof) cfg_line
         if(cfg_line(1:1).eq.'*' .or.  cfg_line(2:2).eq.'*') exit
         if (eof /= 0) exit
         if(cfg_line(3:3).eq.' '.or.cfg_line(4:4).eq.' ') exit
!         if(cfg_line(3:3).eq.' '.or.cfg_line(4:4).eq.' ') THEN
!           print*, "STOP analy1: A"
!           STOP
!         END IF
!         if(cfg_line(1:1).ne.' '.or.cfg_line(9:9).ne.' ') exit
         if(cfg_line(1:1).ne.' '.or.cfg_line(9:9).ne.' ') THEN
           print*, "STOP analy1: B"
           STOP
         END IF
         read (iread_csf,'(a72)') line
         write (1000+iread_csf,'(a55)') cfg_line
         write (1000+iread_csf,'(a72)') line
         if(I_CSF == 0) THEN
            if(line(1:5).ne.'    '.or.line(6:6).eq.' ') exit
            if(line(7:7).eq.' '.or.line(8:8).eq.' ') exit
         else if(I_CSF == 1) THEN
!            if(line(1:1).ne.' '.or.line(2:2).eq.' ') exit
            if(line(1:1).ne.' '.or.line(2:2).eq.' ') THEN
           print*, "STOP analy1: E"
           STOP
         END IF
!            if(line(3:3).eq.' '.or.line(4:4).eq.' ') exit
            if(line(3:3).eq.' '.or.line(4:4).eq.' ') THEN
           print*, "STOP analy1: F"
           STOP
         END IF
         else
            write(*,*) ' error in cfg.inp format I_CSF',I_CSF
            stop
         end if
         ilinecount = ilinecount + 1
      end do
      rewind(1000+iread_csf)
      if (ilinecount == 0) then
         write(*,*) ' error in cfg.inp format'
         stop
      end if
      csfs_LS%nr_of_csfs=ilinecount
!      write(*,*) '   end subroutine analy1'
      end subroutine analy1

!
!***********************************************************************
!                                                                      *
      subroutine cfgo1(I_CSF,cfg_line)
!                                                                      *
!     This subroutine  up the arrays                                   *
!     csfs_LS%csfs, with data from iread_csf                           *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: I_CSF
      character(len=72), intent(in) :: cfg_line
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=72) :: cfg_line2, csf_line
      character(len=8)  :: somech
      character(len=1)  :: somech1
      character         :: charl
      integer :: lval
      integer :: imultiplicity, itermnr
      integer :: icsf, isuba, isubc, j, jj, i, ii
      integer :: inosubc_count
      integer :: shell_n_1, shell_l_1, shell_n_2, shell_l_2
      integer :: shell_n_3, shell_l_3, shell_n_4, shell_l_4
      integer :: shell_n_5, shell_l_5, shell_n_6, shell_l_6
      integer :: shell_n_7, shell_l_7, shell_n_8, shell_l_8
      integer :: shell_n_9, shell_l_9, shell_n_10, shell_l_10
      integer :: shell_n_11, shell_l_11, shell_n_12, shell_l_12
      integer :: iN_T, iL_T, iN_big_T, iJ1, iJ2
      integer, dimension(2) :: iJJ
!-----------------------------------------------
      shell_n_1  = -100
      shell_n_2  = -100
      shell_n_3  = -100
      shell_n_4  = -100
      shell_n_5  = -100
      shell_n_6  = -100
      shell_n_7  = -100
      shell_n_8  = -100
      shell_n_9  = -100
      shell_n_10 = -100
      shell_n_11 = -100
      shell_n_12 = -100
      allocate (Iselect(csfs_LS%nr_of_csfs))
!      write(*,*) 'Specify shells for recoupling (no more than 12)'
!      read(*,'(a56)') cfg_line
!      write(*,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(*,*) cfg_line
!      write(*,*) ""
!      write(iwrite_log,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(iwrite_log,*) cfg_line
!      write(iwrite_log,*)""
!      write(iwrite_log,*)                                 &
!      "---------------------------------------------------------------"
!      write(iwrite_log,*)"List of configurations which",  &
!      " are eliminated from the calculation"
!      write(iwrite_log,*)''
!      write(iwrite_log,*)"   No.           Configurations"
!      write(iwrite_log,*)                                 &
!      "---------------------------------------------------------------"
      read(cfg_line(1:1),'(a1)') somech1
      if (somech1(1:1) == ' ') then
         if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(2:2),'(i1)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         else if(cfg_line(5:5) == ' ' .or. cfg_line(5:5) .eq. ',') then
            read(cfg_line(2:3),'(i2)') shell_n_1
            read(cfg_line(4:4),'(a1)') somech1
            jj=5
         end if
      else
         if(cfg_line(3:3) == ' ' .or. cfg_line(3:3) .eq. ',') then
            read(cfg_line(1:1),'(i1)') shell_n_1
            read(cfg_line(2:2),'(a1)') somech1
            jj=3
         else if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(1:2),'(i2)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         end if
      end if
      shell_l_1 = lval(somech1)
    1 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
!         write(*,*) 'Must be two shells for recoupling: cfgo2'
!         stop
!      end if
       if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
         jj=jj+1
         go to 1
       else
         if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_2
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
             jj=jj+2
         else if(cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
           read(cfg_line(jj:jj+1),'(i2)') shell_n_2
           read(cfg_line(jj+2:jj+2),'(a1)') somech1
           jj=jj+3
         end if
       end if
       shell_l_2 = lval(somech1)
      end if
    2 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 2
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_3
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_3
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_3 = lval(somech1)
      end if
    3 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 3
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_4
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_4
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_4 = lval(somech1)
      end if
    4 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 4
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_5
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_5
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_5 = lval(somech1)
      end if
    5 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 5
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_6
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_6
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_6 = lval(somech1)
      end if
    6 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 6
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_7
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_7
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_7 = lval(somech1)
      end if
    7 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 7
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_8
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_8
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_8 = lval(somech1)
      end if
    8 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 8
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_9
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_9
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_9 = lval(somech1)
      end if
    9 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 9
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_10
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_10
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_10 = lval(somech1)
      end if
   10 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 10
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_11
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_11
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_11 = lval(somech1)
      end if
   11 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 11
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
             read(cfg_line(jj:jj),'(i1)') shell_n_12
             read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_12
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_12 = lval(somech1)
      end if
!
!      read(iread_csf,'(a1)')csf_line
!      read(iread_csf,'(a1)')csf_line
!zr froese aprasyme ka reiskia tie numeriai prie termu...
      do i=1,72
         cfg_line2(i:i)=' '
         csf_line(i:i)=' '
      end do
!write(*,*)'   ', cfgs_LS%nr_of_cfg
      do icsf=1, csfs_LS%nr_of_csfs,1
         inosubc_count=0
         read(1000+iread_csf,'(a56)')cfg_line2
         read(1000+iread_csf,'(a72)')csf_line
!define the number of coupled shells nosubc
         do j=1, 7, 1
            read(cfg_line2(8*j-7:8*j),'(a8)') somech
            if(somech.ne.'      ') then
!write(*,*) '*',somech,'*'
              inosubc_count=inosubc_count+1
            end if
!---- sitas apsisaugojimas tam atvejui,
!jei tai ne cfg.inp, o cfg.out failas ...
            if(somech.eq.'      ') exit
         end do
!GG         csfs_LS%csfs(icsf)%nosubc=inosubc_count
         csfs_LS%csfs(icsf)%nosubc=1
!write(*,*)'inosubc_count=', inosubc_count
         allocate &
             (csfs_LS%csfs(icsf)%subc_cfg(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%subc(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM1(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM2(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iJ(csfs_LS%csfs(icsf)%nosubc))
!GG Cia kaszas keisto!!!         allocate(csfs_LS%csfs(icsf)%iJ(csfs_LS%csfs(icsf)%nosubc))
!uzpildom cfg ...
!read quantum numbers of each subshell
!         do isubc=1,csfs_LS%csfs(icsf)%nosubc,1
         isuba = 0
         do isubc=1,inosubc_count,1
!read quantum numbers N, n, and l of csf configuration
            read(cfg_line2(8*isubc-6:8*isubc-5),'(i2)') iN_T
            read(cfg_line2(8*isubc-4:8*isubc-4),'(a1)') charl
            iL_T= lval(charl)
            read(cfg_line2(8*isubc-2:8*isubc-1),'(i2)') iN_big_T
            if(iN_T .eq. shell_n_1 .and. iL_T .eq. shell_l_1) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_2 .and. iL_T .eq. shell_l_2) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_3 .and. iL_T .eq. shell_l_3) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_4 .and. iL_T .eq. shell_l_4) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_5 .and. iL_T .eq. shell_l_5) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_6 .and. iL_T .eq. shell_l_6) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_7 .and. iL_T .eq. shell_l_7) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_8 .and. iL_T .eq. shell_l_8) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_9 .and. iL_T .eq. shell_l_9) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_10 .and. iL_T .eq. shell_l_10) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_11 .and. iL_T .eq. shell_l_11) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_12 .and. iL_T .eq. shell_l_12) then
               isuba = isuba + 1
            else
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            if (isuba == 1) then
              Iselect(icsf) = 1
            else if (isuba >= 2) then
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            iJJ(isuba) = isubc
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%iN_big = iN_big_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%in = iN_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%il= iL_T
            if(I_CSF == 0) THEN
               read(csf_line(8*isubc-7:8*isubc),'(a8)') somech
               read(somech(5:6),'(I2)')imultiplicity
               read(somech(7:7),'(a1)')charl
               read(somech(8:8),'(I1)')itermnr
            else if(I_CSF == 1) THEN
               read(csf_line(4*isubc-3:4*isubc),'(a8)') somech
               read(somech(1:2),'(I2)')imultiplicity
               read(somech(3:3),'(a1)')charl
               read(somech(4:4),'(I1)')itermnr
            end if
            if(imultiplicity.le.0) then
               write(iwrite_log,'(a22,i4,a4,i1,a5,i2)')     &
               ' error reading csf Nr. ',icsf,' 2*S',isubc, &
               '+1 = ',imultiplicity
               stop 'unrecognized csf Si'
            end if
            if(lval(charl).le.(-1)) then
               write(iwrite_log,'(a22,i4,a2,i1,a5,i2)')       &
               ' error reading csf Nr. ',icsf,' L',isubc,'= ',&
               lval(charl)
               stop 'unrecognized csf Li'
            end if
            if(itermnr.le.(-1)) then
               write(iwrite_log,'(a22,i4,a3,i1,a5,i2)') &
               ' error reading csf Nr. ',icsf,' nr',isubc,'= ',itermnr
               stop 'unrecognized csf nri'
            end if
            csfs_LS%csfs(icsf)%subc(isuba)%iS=imultiplicity-1
            csfs_LS%csfs(icsf)%subc(isuba)%iL=lval(charl)*2
            csfs_LS%csfs(icsf)%subc(isuba)%inr=itermnr
         end do
         if (isuba /= 1) then
           Iselect(icsf) = 0
         end if
         if (Iselect(icsf) == 0) then
            write(iwrite_log,'(I7,A)') icsf, cfg_line2
            cycle
         end if
         csfs_LS%csfs(icsf)%iM1(1)=csfs_LS%csfs(icsf)%subc(1)%iL
         csfs_LS%csfs(icsf)%iM2(1)=csfs_LS%csfs(icsf)%subc(1)%iS
         j = iJJ(2)
         if(I_CSF == 0) THEN
           read(                                                       &
           csf_line(8*inosubc_count+8*(j-1)-7:8*inosubc_count+8*(j-1)),&
           '(a8)') somech
           read(somech(5:6),'(I2)')imultiplicity
           read(somech(7:7),'(a1)')charl
         else if(I_CSF == 1) THEN
           read(                                                       &
           csf_line(4*inosubc_count+4*(j-1)-3:4*inosubc_count+4*(j-1)),&
           '(a8)') somech
           read(somech(1:2),'(I2)')imultiplicity
           read(somech(3:3),'(a1)')charl
         end if
         csfs_LS%csfs(icsf)%iM1(2)=lval(charl)*2
         csfs_LS%csfs(icsf)%iM2(2)=imultiplicity-1
      end do
      write(iwrite_log,*)                                 &
      "---------------------------------------------------------------"
!      write(*,*) '   end subroutine cfgo1'
      end subroutine cfgo1
!
!***********************************************************************
!                                                                      *
      subroutine cfgo2(I_CSF,cfg_line)
!                                                                      *
!     This subroutine  up the arrays                                   *
!     csfs_LS%csfs, with data from iread_csf                           *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: I_CSF
      character(len=72), intent(in) :: cfg_line
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=72) :: cfg_line2, csf_line
      character(len=8)  :: somech
      character(len=1)  :: somech1
      character    :: charl
      integer :: lval
      integer :: imultiplicity, itermnr
      integer :: icsf, isuba, isubc, j, jj, i, ii
      integer :: inosubc_count
      integer :: shell_n_1, shell_l_1, shell_n_2, shell_l_2
      integer :: shell_n_3, shell_l_3, shell_n_4, shell_l_4
      integer :: shell_n_5, shell_l_5, shell_n_6, shell_l_6
      integer :: shell_n_7, shell_l_7, shell_n_8, shell_l_8
      integer :: shell_n_9, shell_l_9, shell_n_10, shell_l_10
      integer :: shell_n_11, shell_l_11, shell_n_12, shell_l_12
      integer :: iN_T, iL_T, iN_big_T, iJ1, iJ2
      integer, dimension(2) :: iJJ
!-----------------------------------------------
      shell_n_1  = -100
      shell_n_2  = -100
      shell_n_3  = -100
      shell_n_4  = -100
      shell_n_5  = -100
      shell_n_6  = -100
      shell_n_7  = -100
      shell_n_8  = -100
      shell_n_9  = -100
      shell_n_10 = -100
      shell_n_11 = -100
      shell_n_12 = -100
      allocate (Iselect(csfs_LS%nr_of_csfs))
!      write(*,*) 'Specify shells for recoupling (no more than 12)'
!      read(*,'(a56)') cfg_line
!      write(*,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(*,*) cfg_line
!      write(*,*) ""
!      write(iwrite_log,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(iwrite_log,*) cfg_line
!      write(iwrite_log,*)""
!      write(iwrite_log,*)                                 &
!      "---------------------------------------------------------------"
!      write(iwrite_log,*)"List of configurations which",  &
!      " are eliminated from the calculation"
!      write(iwrite_log,*)''
!      write(iwrite_log,*)"   No.           Configurations"
!      write(iwrite_log,*)                                 &
!      "------  -------------------------------------------------------"
      read(cfg_line(1:1),'(a1)') somech1
      if (somech1(1:1) == ' ') then
         if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(2:2),'(i1)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         else if(cfg_line(5:5) == ' ' .or. cfg_line(5:5) .eq. ',') then
            read(cfg_line(2:3),'(i2)') shell_n_1
            read(cfg_line(4:4),'(a1)') somech1
            jj=5
         end if
      else
         if(cfg_line(3:3) == ' ' .or. cfg_line(3:3) .eq. ',') then
            read(cfg_line(1:1),'(i1)') shell_n_1
            read(cfg_line(2:2),'(a1)') somech1
            jj=3
         else if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(1:2),'(i2)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         end if
      end if
      shell_l_1 = lval(somech1)
    1 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj >= 57) then
         write(*,*) 'Must be two shells for recoupling: cfgo1'
         stop
      end if
      if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
        jj=jj+1
        go to 1
      else
        if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
           read(cfg_line(jj:jj),'(i1)') shell_n_2
           read(cfg_line(jj+1:jj+1),'(a1)') somech1
           jj=jj+2
        else if(cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
           read(cfg_line(jj:jj+1),'(i2)') shell_n_2
           read(cfg_line(jj+2:jj+2),'(a1)') somech1
           jj=jj+3
        end if
      end if
      shell_l_2 = lval(somech1)
    2 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 2
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_3
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_3
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_3 = lval(somech1)
      end if
    3 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 3
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_4
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_4
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_4 = lval(somech1)
      end if
    4 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 4
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_5
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_5
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_5 = lval(somech1)
      end if
    5 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 5
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_6
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_6
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_6 = lval(somech1)
      end if
    6 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 6
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_7
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_7
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_7 = lval(somech1)
      end if
    7 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 7
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_8
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_8
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_8 = lval(somech1)
      end if
    8 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 8
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_9
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_9
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_9 = lval(somech1)
      end if
    9 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 9
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_10
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_10
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_10 = lval(somech1)
      end if
   10 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 10
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_11
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_11
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_11 = lval(somech1)
      end if
   11 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 11
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_12
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_12
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_12 = lval(somech1)
      end if
!
!      read(iread_csf,'(a1)')csf_line
!      read(iread_csf,'(a1)')csf_line
!zr froese aprasyme ka reiskia tie numeriai prie termu...
      do i=1,72
         cfg_line2(i:i)=' '
         csf_line(i:i)=' '
      end do
!      write(*,*)'cfgs_LS%nr_of_cfg   ', cfgs_LS%nr_of_cfg
      do icsf=1, csfs_LS%nr_of_csfs,1
         inosubc_count=0
         read(1000+iread_csf,'(a56)')cfg_line2
         read(1000+iread_csf,'(a72)')csf_line
!define the number of coupled shells nosubc
         do j=1, 7, 1
            read(cfg_line2(8*j-7:8*j),'(a8)') somech
            if(somech.ne.'      ') then
!write(*,*) '*',somech,'*'
              inosubc_count=inosubc_count+1
            end if
!---- sitas apsisaugojimas tam atvejui,
!jei tai ne cfg.inp, o cfg.out failas ...
            if(somech.eq.'      ') exit
         end do
!GG         csfs_LS%csfs(icsf)%nosubc=inosubc_count
         csfs_LS%csfs(icsf)%nosubc=2
!write(*,*)'inosubc_count=', inosubc_count
         allocate &
             (csfs_LS%csfs(icsf)%subc_cfg(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%subc(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM1(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM2(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iJ(csfs_LS%csfs(icsf)%nosubc))
!GG Cia kaszas keisto!!!         allocate(csfs_LS%csfs(icsf)%iJ(csfs_LS%csfs(icsf)%nosubc))
!uzpildom cfg ...
!read quantum numbers of each subshell
!         do isubc=1,csfs_LS%csfs(icsf)%nosubc,1
         isuba = 0
         do isubc=1,inosubc_count,1
!read quantum numbers N, n, and l of csf configuration
            read(cfg_line2(8*isubc-6:8*isubc-5),'(i2)') iN_T
            read(cfg_line2(8*isubc-4:8*isubc-4),'(a1)') charl
            iL_T= lval(charl)
            read(cfg_line2(8*isubc-2:8*isubc-1),'(i2)') iN_big_T
            if(iN_T .eq. shell_n_1 .and. iL_T .eq. shell_l_1) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_2 .and. iL_T .eq. shell_l_2) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_3 .and. iL_T .eq. shell_l_3) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_4 .and. iL_T .eq. shell_l_4) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_5 .and. iL_T .eq. shell_l_5) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_6 .and. iL_T .eq. shell_l_6) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_7 .and. iL_T .eq. shell_l_7) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_8 .and. iL_T .eq. shell_l_8) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_9 .and. iL_T .eq. shell_l_9) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_10 .and. iL_T .eq. shell_l_10) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_11 .and. iL_T .eq. shell_l_11) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_12 .and. iL_T .eq. shell_l_12) then
               isuba = isuba + 1
            else
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            if (isuba == 1) then
              Iselect(icsf) = 1
            else if (isuba >= 3) then
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            iJJ(isuba) = isubc
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%iN_big = iN_big_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%in = iN_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%il= iL_T
            if(I_CSF == 0) THEN
               read(csf_line(8*isubc-7:8*isubc),'(a8)') somech
               read(somech(5:6),'(I2)')imultiplicity
               read(somech(7:7),'(a1)')charl
               read(somech(8:8),'(I1)')itermnr
            else if(I_CSF == 1) THEN
               read(csf_line(4*isubc-3:4*isubc),'(a8)') somech
               read(somech(1:2),'(I2)')imultiplicity
               read(somech(3:3),'(a1)')charl
               read(somech(4:4),'(I1)')itermnr
            end if
            if(imultiplicity.le.0) then
               write(iwrite_log,'(a22,i4,a4,i1,a5,i2)')     &
               ' error reading csf Nr. ',icsf,' 2*S',isubc, &
               '+1 = ',imultiplicity
               stop 'unrecognized csf Si'
            end if
            if(lval(charl).le.(-1)) then
               write(iwrite_log,'(a22,i4,a2,i1,a5,i2)')       &
               ' error reading csf Nr. ',icsf,' L',isubc,'= ',&
               lval(charl)
               stop 'unrecognized csf Li'
            end if
            if(itermnr.le.(-1)) then
               write(iwrite_log,'(a22,i4,a3,i1,a5,i2)') &
               ' error reading csf Nr. ',icsf,' nr',isubc,'= ',itermnr
               stop 'unrecognized csf nri'
            end if
            csfs_LS%csfs(icsf)%subc(isuba)%iS=imultiplicity-1
            csfs_LS%csfs(icsf)%subc(isuba)%iL=lval(charl)*2
            csfs_LS%csfs(icsf)%subc(isuba)%inr=itermnr
         end do
         if (isuba == 1) then
           Iselect(icsf) = 0
         end if
         if (Iselect(icsf) == 0) then
            write(iwrite_log,'(I7,A)') icsf, cfg_line2
            cycle
         end if
         csfs_LS%csfs(icsf)%iM1(1)=csfs_LS%csfs(icsf)%subc(1)%iL
         csfs_LS%csfs(icsf)%iM2(1)=csfs_LS%csfs(icsf)%subc(1)%iS
         j = iJJ(2)
         if(I_CSF == 0) THEN
           read(                                                       &
           csf_line(8*inosubc_count+8*(j-1)-7:8*inosubc_count+8*(j-1)),&
           '(a8)') somech
           read(somech(5:6),'(I2)')imultiplicity
           read(somech(7:7),'(a1)')charl
         else if(I_CSF == 1) THEN
           read(                                                       &
           csf_line(4*inosubc_count+4*(j-1)-3:4*inosubc_count+4*(j-1)),&
           '(a8)') somech
           read(somech(1:2),'(I2)')imultiplicity
           read(somech(3:3),'(a1)')charl
         end if
         csfs_LS%csfs(icsf)%iM1(2)=lval(charl)*2
         csfs_LS%csfs(icsf)%iM2(2)=imultiplicity-1
      end do
      write(iwrite_log,*)                                 &
      "---------------------------------------------------------------"
!      write(*,*) '   end subroutine cfgo2'
      end subroutine cfgo2
!
!***********************************************************************
!                                                                      *
      subroutine cfgo3(I_CSF,cfg_line)
!                                                                      *
!     This subroutine  up the arrays                                   *
!     csfs_LS%csfs, with data from iread_csf                           *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: I_CSF
      character(len=72), intent(in) :: cfg_line
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=72) :: cfg_line2, csf_line
      character(len=8)  :: somech
      character(len=1)  :: somech1
      character    :: charl
      integer :: lval
      integer :: imultiplicity, itermnr
      integer :: icsf, isuba, isubc, j, jj, i, ii
      integer :: inosubc_count
      integer :: shell_n_1, shell_l_1, shell_n_2, shell_l_2
      integer :: shell_n_3, shell_l_3, shell_n_4, shell_l_4
      integer :: shell_n_5, shell_l_5, shell_n_6, shell_l_6
      integer :: shell_n_7, shell_l_7, shell_n_8, shell_l_8
      integer :: shell_n_9, shell_l_9, shell_n_10, shell_l_10
      integer :: shell_n_11, shell_l_11, shell_n_12, shell_l_12
      integer :: iN_T, iL_T, iN_big_T, iJ1, iJ2
      integer, dimension(3) :: iJJ
!-----------------------------------------------
      shell_n_1  = -100
      shell_n_2  = -100
      shell_n_3  = -100
      shell_n_4  = -100
      shell_n_5  = -100
      shell_n_6  = -100
      shell_n_7  = -100
      shell_n_8  = -100
      shell_n_9  = -100
      shell_n_10 = -100
      shell_n_11 = -100
      shell_n_12 = -100
      allocate (Iselect(csfs_LS%nr_of_csfs))
!      write(*,*) 'Specify shells for recoupling (no more than 12)'
!      read(*,'(a56)') cfg_line
!      write(*,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(*,*) cfg_line
!      write(*,*) ""
!      write(iwrite_log,*)""
!      write(iwrite_log,*)                                 &
!                      "List of the shells included in the recoupling:"
!      write(iwrite_log,*) cfg_line
!      write(iwrite_log,*)""
!      write(iwrite_log,*)                                 &
!      "---------------------------------------------------------------"
!      write(iwrite_log,*)"List of configurations which",  &
!      " are eliminated from the calculation"
!      write(iwrite_log,*)''
!      write(iwrite_log,*)"   No.           Configurations"
!      write(iwrite_log,*)                                 &
!      "------  -------------------------------------------------------"
      read(cfg_line(1:1),'(a1)') somech1
      if (somech1(1:1) == ' ') then
         if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(2:2),'(i1)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         else if(cfg_line(5:5) == ' ' .or. cfg_line(5:5) .eq. ',') then
            read(cfg_line(2:3),'(i2)') shell_n_1
            read(cfg_line(4:4),'(a1)') somech1
            jj=5
         end if
      else
         if(cfg_line(3:3) == ' ' .or. cfg_line(3:3) .eq. ',') then
            read(cfg_line(1:1),'(i1)') shell_n_1
            read(cfg_line(2:2),'(a1)') somech1
            jj=3
         else if(cfg_line(4:4) == ' ' .or. cfg_line(4:4) .eq. ',') then
            read(cfg_line(1:2),'(i2)') shell_n_1
            read(cfg_line(3:3),'(a1)') somech1
            jj=4
         end if
      end if
      shell_l_1 = lval(somech1)
    1 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj >= 57) then
         write(*,*) 'Must be two shells for recoupling: cfgo3'
         stop
      end if
      if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
        jj=jj+1
        go to 1
      else
        if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
           read(cfg_line(jj:jj),'(i1)') shell_n_2
           read(cfg_line(jj+1:jj+1),'(a1)') somech1
           jj=jj+2
        else if(cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
           read(cfg_line(jj:jj+1),'(i2)') shell_n_2
           read(cfg_line(jj+2:jj+2),'(a1)') somech1
           jj=jj+3
        end if
      end if
      shell_l_2 = lval(somech1)
    2 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 2
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_3
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_3
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_3 = lval(somech1)
      end if
    3 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 3
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_4
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_4
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_4 = lval(somech1)
      end if
    4 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 4
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_5
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_5
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_5 = lval(somech1)
      end if
    5 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 5
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_6
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_6
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_6 = lval(somech1)
      end if
    6 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 6
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_7
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_7
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_7 = lval(somech1)
      end if
    7 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 7
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_8
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_8
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_8 = lval(somech1)
      end if
    8 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 8
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_9
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_9
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_9 = lval(somech1)
      end if
    9 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 9
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_10
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_10
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_10 = lval(somech1)
      end if
   10 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 10
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_11
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_11
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_11 = lval(somech1)
      end if
   11 continue
      read(cfg_line(jj:jj),'(a1)') somech1
      if (jj <= 57) then
         if (somech1 .eq. ' ' .or. somech1 .eq. ',') then
            jj=jj+1
            go to 11
         else
           if(cfg_line(jj+2:jj+2)==' '.or.cfg_line(jj+2:jj+2).eq.',') then
              read(cfg_line(jj:jj),'(i1)') shell_n_12
              read(cfg_line(jj+1:jj+1),'(a1)') somech1
              jj=jj+2
            else if                                                   &
              (cfg_line(jj+3:jj+3)==' '.or.cfg_line(jj+3:jj+3).eq.',')&
                                                                then
              read(cfg_line(jj:jj+1),'(i2)') shell_n_12
              read(cfg_line(jj+2:jj+2),'(a1)') somech1
              jj=jj+3
            end if
         end if
         shell_l_12 = lval(somech1)
      end if
!
!      read(iread_csf,'(a56)')csf_line
!      read(iread_csf,'(a56)')csf_line
!zr froese aprasyme ka reiskia tie numeriai prie termu...
      do i=1,72
         cfg_line2(i:i)=' '
         csf_line(i:i)=' '
      end do
!write(*,*)'   ', cfgs_LS%nr_of_cfg
      do icsf=1, csfs_LS%nr_of_csfs,1
         read(1000+iread_csf,'(a56)')cfg_line2
         read(1000+iread_csf,'(a72)')csf_line
!define the number of coupled shells nosubc
         inosubc_count=0
         do j=1, 7, 1
            read(cfg_line2(8*j-7:8*j),'(a8)') somech
            if(somech.ne.'      ') then
!write(*,*) '*',somech,'*'
              inosubc_count=inosubc_count+1
            end if
!---- sitas apsisaugojimas tam atvejui,
!jei tai ne cfg.inp, o cfg.out failas ...
            if(somech.eq.'      ') exit
         end do
!GG         csfs_LS%csfs(icsf)%nosubc=inosubc_count
         csfs_LS%csfs(icsf)%nosubc=3
!write(*,*)'inosubc_count=', inosubc_count
         allocate &
             (csfs_LS%csfs(icsf)%subc_cfg(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%subc(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM1(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iM2(csfs_LS%csfs(icsf)%nosubc))
         allocate(csfs_LS%csfs(icsf)%iJ(csfs_LS%csfs(icsf)%nosubc))
!uzpildom cfg ...
!read quantum numbers of each subshell
!         do isubc=1,csfs_LS%csfs(icsf)%nosubc,1
         isuba = 0
         do isubc=1,inosubc_count,1
!read quantum numbers N, n, and l of csf configuration
            read(cfg_line2(8*isubc-6:8*isubc-5),'(i2)') iN_T
            read(cfg_line2(8*isubc-4:8*isubc-4),'(a1)') charl
            iL_T= lval(charl)
            read(cfg_line2(8*isubc-2:8*isubc-1),'(i2)') iN_big_T
            if(iN_T .eq. shell_n_1 .and. iL_T .eq. shell_l_1) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_2 .and. iL_T .eq. shell_l_2) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_3 .and. iL_T .eq. shell_l_3) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_4 .and. iL_T .eq. shell_l_4) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_5 .and. iL_T .eq. shell_l_5) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_6 .and. iL_T .eq. shell_l_6) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_7 .and. iL_T .eq. shell_l_7) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_8 .and. iL_T .eq. shell_l_8) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_9 .and. iL_T .eq. shell_l_9) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_10 .and. iL_T .eq. shell_l_10) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_11 .and. iL_T .eq. shell_l_11) then
               isuba = isuba + 1
            else if(iN_T .eq. shell_n_12 .and. iL_T .eq. shell_l_12) then
               isuba = isuba + 1
            else
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            if (isuba == 1) then
              Iselect(icsf) = 1
            else if (isuba >= 4) then
               if(iN_big_T /= 4*iL_T+2) then
                  Iselect(icsf) = 0
                  exit
               else
                  cycle
               end if
            end if
            iJJ(isuba) = isubc
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%iN_big = iN_big_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%in = iN_T
            csfs_LS%csfs(icsf)%subc_cfg(isuba)%il= iL_T
            if(I_CSF == 0) THEN
               read(csf_line(8*isubc-7:8*isubc),'(a8)') somech
               read(somech(5:6),'(I2)')imultiplicity
               read(somech(7:7),'(a1)')charl
               read(somech(8:8),'(I1)')itermnr
            else if(I_CSF == 1) THEN
               read(csf_line(4*isubc-3:4*isubc),'(a8)') somech
               read(somech(1:2),'(I2)')imultiplicity
               read(somech(3:3),'(a1)')charl
               read(somech(4:4),'(I1)')itermnr
            end if
            if(imultiplicity.le.0) then
               write(iwrite_log,'(a22,i4,a4,i1,a5,i2)')     &
               ' error reading csf Nr. ',icsf,' 2*S',isubc, &
               '+1 = ',imultiplicity
               stop 'unrecognized csf Si'
            end if
            if(lval(charl).le.(-1)) then
               write(iwrite_log,'(a22,i4,a2,i1,a5,i2)')       &
               ' error reading csf Nr. ',icsf,' L',isubc,'= ',&
               lval(charl)
               stop 'unrecognized csf Li'
            end if
            if(itermnr.le.(-1)) then
               write(iwrite_log,'(a22,i4,a3,i1,a5,i2)') &
               ' error reading csf Nr. ',icsf,' nr',isubc,'= ',itermnr
               stop 'unrecognized csf nri'
            end if
            csfs_LS%csfs(icsf)%subc(isuba)%iS=imultiplicity-1
            csfs_LS%csfs(icsf)%subc(isuba)%iL=lval(charl)*2
            csfs_LS%csfs(icsf)%subc(isuba)%inr=itermnr
         end do
         if (isuba == 1 .or. isuba == 2) then
           Iselect(icsf) = 0
         end if
         if (Iselect(icsf) == 0) then
            write(iwrite_log,'(I7,A)') icsf, cfg_line2
            cycle
         end if
         do i =1,3
            j = iJJ(i)
            if(j == 1) then
               csfs_LS%csfs(icsf)%iM1(1)=csfs_LS%csfs(icsf)%subc(1)%iL
               csfs_LS%csfs(icsf)%iM2(1)=csfs_LS%csfs(icsf)%subc(1)%iS
            else
               if(I_CSF == 0) THEN
                  read(csf_line(                                       &
                    8*inosubc_count+8*(j-1)-7:8*inosubc_count+8*(j-1)),&
                    '(a8)') somech
                  read(somech(5:6),'(I2)')imultiplicity
                  read(somech(7:7),'(a1)')charl
               else if(I_CSF == 1) THEN
                  read(csf_line(                                       &
                    4*inosubc_count+4*(j-1)-3:4*inosubc_count+4*(j-1)),&
                    '(a8)') somech
                  read(somech(1:2),'(I2)')imultiplicity
                  read(somech(3:3),'(a1)')charl
               end if
               csfs_LS%csfs(icsf)%iM1(i)=lval(charl)*2
               csfs_LS%csfs(icsf)%iM2(i)=imultiplicity-1
            end if
         end do
      end do
      write(iwrite_log,*)                                 &
      "---------------------------------------------------------------"
!      write(*,*) '   end subroutine cfgo3'
      end subroutine cfgo3
!
!***********************************************************************
!                                                                      *
      subroutine analy_j(icase,I_CSF, icoupling_nr_in_expansions)
!                                                                      *
!     This subroutine analyzes .j file and defines number of           *
!     J (nr_of_j) and number of states (states%nr_of_states)           *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent (in) :: I_CSF ,icoupling_nr_in_expansions,icase
!      integer, intent(out) :: nr_of_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=2) :: atom
      real        :: z
      integer     :: n_of_j,ncfg
      integer     :: numb_of_rows, J, inumber, eof
!-----------------------------------------------
!      write(*,*) '   subroutine analy_j'
      if(icase == 1 .and. I_CSF == 0) THEN
         read(iread_coefs,'(2x,a2,10x,f5.1,7x,i2,11x,i2)') &
                  atom, z, n_of_j, ncfg
      else if(icase == 1 .and. I_CSF == 1) THEN
         read(iread_coefs,'(2x,a2,10x,f5.1,7x,i4,9x,i4)') &
                  atom, z, n_of_j, ncfg
      end if
!      if(ncfg == 0) then
         ncfg = csfs_LS%nr_of_csfs
!      else if(ncfg.ne.csfs_LS%nr_of_csfs) then
!         write(*,*) 'ncfg from j-file does not mach number of csf ', &
!         'defined from csf-file'
!         write(*,*)  'ncfg,csfs_LS%nr_of_csfs',ncfg,csfs_LS%nr_of_csfs
!         stop
!      end if
      if((ncfg/7)*7 == ncfg) then
         numb_of_rows=int(ncfg/7)
      else
         numb_of_rows=int(ncfg/7)+1
      end if
!      nr_of_j=0
!      states%nr_of_states=0
!      do
         read(iread_coefs,*,IOSTAT=eof)
!         if (eof /= 0) exit
         if (eof /= 0) THEN
           print*, "STOP analy_j: A"
           STOP
         END IF
         read(iread_coefs,*,iostat=eof)
!         if (eof /= 0) exit
         if (eof /= 0) THEN
           print*, "STOP analy_j: B"
           STOP
         END IF
         read(iread_coefs,'(2X,I1,7X,I2,11X,I3)',iostat=eof) &
                                             isome, J, inumber
         write(1000+iread_coefs,'(2X,I1,7X,I2,11X,I3)',iostat=eof) &
                                             isome, J, inumber
!         if (eof /= 0) exit
         if (eof /= 0) THEN
           print*, "STOP analy_j: C"
           STOP
         END IF
!         if(isome.ne.2) exit
         if(isome.ne.2) THEN
           print*, "STOP analy_j: D"
           STOP
         END IF
!         nr_of_j = nr_of_j +1
!J_structure%nr_of_j = J_structure%nr_of_j +1
!         states%nr_of_states=states%nr_of_states + inumber
         states%nr_of_states= inumber
!         do i1=1,inumber
!            read(iread_coefs,*)
!            read(iread_coefs,*)
!            do i2=1,numb_of_rows
!               read(iread_coefs,*)
!            end do
!         end do
!      end do
      all_expansions%coupling_expansions(icoupling_nr_in_expansions)% &
                     nr_of_expansions                                 &
      = states%nr_of_states
      write(iwrite_log,'(1x,a19,i4)') 'Number of states = ', &
      states%nr_of_states
      write(*,'(1x,a19,i4)') 'Number of states = ', &
      states%nr_of_states
      rewind(1000+iread_coefs)
!      write(*,*) '   end subroutine analy_j'
      end subroutine analy_j
!
!***********************************************************************
!                                                                      *
      subroutine readcoeffs(icase,I_CSF,nr_of_j,icoupling_nr_in_expansions)
!                                                                      *
!     This subroutine reads the quantum numbers of the states          *
!     (energies, J's) and expansion coefficients from the expansions   *
!     file *.j and stores it in the "states" and  "all_expansions"     *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: I_CSF,nr_of_j, icoupling_nr_in_expansions, &
       icase
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=2)  :: atom
      character(len=77) :: line
      real(kind=dp):: z, Calmixmax, mixmax
      integer      :: n_of_j,ncfg, numb_of_rows, numb_of_rest
      integer      :: J, inumber, isome, istate_count
      real(kind=dp):: energy, sum_of_state
      integer      :: i0, i1, i2, i2t7, i3
      real(kind=dp),dimension(:),pointer::temp_coeffs
      integer      :: icsf, iexpansion
!-----------------------------------------------
      istate_count=0
!      write(*,*) '   subroutine readcoeffs'
!      IF(icase == 1) THEN
!         rewind(iread_coefs)
!      if(I_CSF == 0) THEN
!         read(iread_coefs,'(2x,a2,10x,f5.1,7x,i2,11x,i2)') &
!                  atom, z, n_of_j, ncfg
!      else if(I_CSF == 1) THEN
!         read(iread_coefs,'(2x,a2,10x,f5.1,7x,i4,9x,i4)')  &
!                  atom, z, n_of_j, ncfg
!      end if
!      end if
!      if(ncfg == 0) then
         ncfg = csfs_LS%nr_of_csfs
!      else if(ncfg.ne.csfs_LS%nr_of_csfs) then
!         write(*,*)                                                  &
!         'ncfg from j-file does not mach number of csf defined from',&
!         ' csf-file'
!         write(*,*) 'ncfg,csfs_LS%nr_of_csfs',ncfg,csfs_LS%nr_of_csfs
!         stop
!      end if
!      if(ncfg.eq.0) stop 'STOP at subroutine readcoeffs: ncfg.eq.0'
      numb_of_rows=int(ncfg/7)
      numb_of_rest=ncfg-numb_of_rows*7
      allocate(temp_coeffs(ncfg))
!      do i0=1, nr_of_j
!read blank lines
!         read(iread_coefs,*)
!         read(iread_coefs,*)
!read J, number of eigenvalues for given J matrix
         read(1000+iread_coefs,'(2X,I1,7X,I2,11X,I3)') isome, J, inumber
!read eigenvalues and their sets of coefficients
!for given J
         do i1=1,inumber
            mixmax = 0.0
            Calmixmax = 0.0
            read(iread_coefs,*)
            read(iread_coefs,'(6x,f16.9)') energy
            istate_count=istate_count+1
            if(istate_count.gt.states%nr_of_states) &
            stop 'istate_count.gt.states%nr_of_states'
            states%states(istate_count)%energy=energy
            states%states(istate_count)%J=J
            all_expansions%                                           &
                coupling_expansions(icoupling_nr_in_expansions)%      &
                expansions(istate_count)%nr_of_state                  &
            = istate_count
            do i2=1,numb_of_rows
               if(I_CSF == 0) THEN
                  read(iread_coefs, '(a70)') line
               else if(I_CSF == 1) THEN
                  read(iread_coefs, '(a77)') line
               end if
               i2t7=(i2-1)*7
               do i3=1,7
                  if((i2t7+i3).gt.ncfg) then
                     write(*,*) 'STOP at subroutine readcoeffs: ', &
                                '(i2t7+i3).gt.ncfg'
                     stop
                  end if
                     if(I_CSF == 0) THEN
!GG              read(line(i3*10-8:i3*10),'(f9.7)') temp_coeffs(i2t7+i3)
                     read(line(i3*10-9:i3*10),'(f10.7)') &
                                                    temp_coeffs(i2t7+i3)
                  else if(I_CSF == 1) THEN
                     read(line(i3*11-10:i3*11),'(f11.8)') &
                                                    temp_coeffs(i2t7+i3)
                  end if
               end do
            end do
!read last coefs of given eigenvalue
            if(numb_of_rest /= 0) then
              if(I_CSF == 0) THEN
                 read(iread_coefs, '(a70)') line
              else if(I_CSF == 1) THEN
                 read(iread_coefs, '(a77)') line
              end if
              i2=numb_of_rows*7
              do i3=1, numb_of_rest
                 if((i2+i3).gt.ncfg) then
                    write(*,*) 'STOP at subroutine readcoeffs: ', &
                               '(i2+i3).gt.ncfg'
                    stop
                 end if
                 if(I_CSF == 0) THEN
!GG                read(line(i3*10-8:i3*10),'(f9.7)') temp_coeffs(i2+i3)
                  read(line(i3*10-9:i3*10),'(f10.7)') temp_coeffs(i2+i3)
                 else if(I_CSF == 1) THEN
                  read(line(i3*11-10:i3*11),'(f11.8)')temp_coeffs(i2+i3)
                 end if
              end do
            end if
!
!find out the size of expansion (nr of csfs with coeffs
! greater than dp_coeff_precision)
            iexpansion=0
            do icsf=1, ncfg
               if(dabs(temp_coeffs(icsf)).gt.dp_coeff_precision) then
                  if (Iselect(icsf) == 1) then
                     iexpansion = iexpansion + 1
                     if(dabs(temp_coeffs(icsf)) > Calmixmax)     &
                              Calmixmax=dabs(temp_coeffs(icsf))
                  else
                     if(dabs(temp_coeffs(icsf)) > mixmax)        &
                              mixmax=dabs(temp_coeffs(icsf))
                  end if
               end if
            end do
            all_expansions%                                           &
                coupling_expansions(icoupling_nr_in_expansions)%      &
                expansions(istate_count)%size                         &
            = iexpansion
            allocate(all_expansions%                                  &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     expansions(istate_count)%coeffs(all_expansions%  &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     expansions(istate_count)%size))
            allocate(all_expansions%                                  &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     expansions(istate_count)%csfs(all_expansions%    &
                     coupling_expansions(icoupling_nr_in_expansions)% &
                     expansions(istate_count)%size))
            iexpansion=0
            sum_of_state = 0.
            do icsf=1, ncfg, 1
               if(dabs(temp_coeffs(icsf)).gt.dp_coeff_precision) then
                  if (Iselect(icsf) == 1) then
                     iexpansion = iexpansion + 1
                     if(all_expansions%                                &
                       coupling_expansions(icoupling_nr_in_expansions)%&
                          expansions(istate_count)%size.lt.iexpansion) &
                                                                    then
                       write(*,*) 'STOP at subroutine readcoeffs: ',   &
              'all_expansions%coupling_expansions()%expansions()%size',&
                       '.lt.iexpansion'
                       stop
                     end if
                     all_expansions%                                   &
                       coupling_expansions(icoupling_nr_in_expansions)%&
                       expansions(istate_count)%coeffs(iexpansion)     &
                     = temp_coeffs(icsf)
!unsing this I am SURE that in "b=a" "b" is "empty" csf ...
                     all_expansions%                                   &
                       coupling_expansions(icoupling_nr_in_expansions)%&
                       expansions(istate_count)%csfs(iexpansion)       &
                     = csfs_LS%csfs(icsf)
!to check the condition: sum of all coeffs^2 = 1
                     sum_of_state = sum_of_state +                     &
                     all_expansions%                                   &
                       coupling_expansions(icoupling_nr_in_expansions)%&
                       expansions(istate_count)%coeffs(iexpansion)*    &
                     all_expansions%                                   &
                       coupling_expansions(icoupling_nr_in_expansions)%&
                       expansions(istate_count)%coeffs(iexpansion)
                  end if
               end if
            end do
            if(i1 ==1)write(iwrite_log,*)                              &
      '---------------------------------------------------------------------'
            if(i1 ==1) write(iwrite_log,*)                             &
      'List of maximum values of expansion coefficients and  Summation rules'
            if(i1 ==1)write(iwrite_log,*)                              &
      '---------------------------------------------------------------------'
            if(i1 ==1) write(iwrite_log,*)                             &
         'State   Maximum value   Maximum value   Expansion   Summation'
            if(i1 ==1) write(iwrite_log,*)                             &
         'No.     included in     removed from    size        rules'
            if(i1 ==1) write(iwrite_log,*)                             &
         '        calculations    calculations'
            if(i1 ==1)write(iwrite_log,*)                              &
      '---------------------------------------------------------------------'
            if (iexpansion < ncfg) then
               write(iwrite_log,'(I6,2X,f11.8,5X,f11.8,I12,3X,f20.17)')&
               istate_count,Calmixmax,mixmax,ncfg,sum_of_state
            else

            write(iwrite_log,'(I6,29X,I12,3X,f20.17)')                 &
                                         istate_count,ncfg,sum_of_state
            end if
         end do
!      end do
      if(associated(temp_coeffs)) deallocate(temp_coeffs)
      write(iwrite_log,*)                                              &
      '---------------------------------------------------------------------'
      write(iwrite_log,*)""
!      write(*,*) '   end subroutine readcoeffs'
      end subroutine readcoeffs
!
!***********************************************************************
!                                                                      *
      subroutine print_csf
!                                                                      *
!     This sbroutine prints the data to the file                       *
!     (was used for the test)                                          *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer     :: icsf, isubc, jsubc,jjsubc, izero
      integer     :: iN_big, in, il, i
      integer     :: inr, Li, Si, L_i, S_i
      character(len=1) :: CVAL
!-----------------------------------------------
      write(iwrite_log,*)                                 &
      "--------------------------------------------------"
      write(iwrite_log,*) &
      "List of CSF which are included in the calculation:"
      i = 0
      do icsf=1,csfs_LS%nr_of_csfs
         if (Iselect(icsf) == 0) then
            cycle
         end if
         i = i + 1
         write(iwrite_log,*)''
         write(iwrite_log,'(3x,a8,i7,6x,a9,i7)')       &
               'CSFs No.', icsf,                        &
               'Shell No.', csfs_LS%csfs(icsf)%nosubc
         write(iwrite_log,'(6x,7(i2,a1,"(",i1,")"))')                 &
                         (csfs_LS%csfs(icsf)%subc_cfg(jsubc)%in,       &
                         CVAL(1,csfs_LS%csfs(icsf)%subc_cfg(jsubc)%il),&
                         csfs_LS%csfs(icsf)%subc_cfg(jsubc)%iN_big,    &
                         jsubc=1,csfs_LS%csfs(icsf)%nosubc)
         if(csfs_LS%csfs(icsf)%nosubc == 1) then
            write(iwrite_log,'(9x,7(1x,i2,a1,i1,1x))')                &
                                 csfs_LS%csfs(icsf)%subc(1)%iS, &
                                 CVAL(2,csfs_LS%csfs(icsf)%subc(1)%iL),&
                                 csfs_LS%csfs(icsf)%subc(1)%inr
         else if (csfs_LS%csfs(icsf)%nosubc == 2) then
            write(iwrite_log,'(9x,7(1x,i2,a1,i1,1x))')                &
                             (csfs_LS%csfs(icsf)%subc(jsubc)%iS+1,     &
                             CVAL(2,csfs_LS%csfs(icsf)%subc(jsubc)%iL),&
                             csfs_LS%csfs(icsf)%subc(jsubc)%inr,       &
                             jsubc=1,csfs_LS%csfs(icsf)%nosubc),       &
                             csfs_LS%csfs(icsf)%iM2(2)+1,              &
                             CVAL(2,csfs_LS%csfs(icsf)%iM1(2))
         else
            izero = 0
            write(iwrite_log,'(9x,7(1x,i2,a1,i1,1x))')                &
                            (csfs_LS%csfs(icsf)%subc(jsubc)%iS+1,      &
                            CVAL(2,csfs_LS%csfs(icsf)%subc(jsubc)%iL), &
                            csfs_LS%csfs(icsf)%subc(jsubc)%inr,        &
                            jsubc=1,csfs_LS%csfs(icsf)%nosubc),        &
                            (csfs_LS%csfs(icsf)%iM2(jjsubc)+1,         &
                            CVAL(2,csfs_LS%csfs(icsf)%iM1(jjsubc)),    &
                            izero,                                     &
                            jjsubc=2,csfs_LS%csfs(icsf)%nosubc)
         end if
         do isubc=1,csfs_LS%csfs(icsf)%nosubc,1
            iN_big= csfs_LS%csfs(icsf)%subc_cfg(isubc)%iN_big
            il = csfs_LS%csfs(icsf)%subc_cfg(isubc)%il
            in = csfs_LS%csfs(icsf)%subc_cfg(isubc)%in
         end do
         do isubc=1,csfs_LS%csfs(icsf)%nosubc,1
            inr= csfs_LS%csfs(icsf)%subc(isubc)%inr
            Li = csfs_LS%csfs(icsf)%subc(isubc)%iL
            Si = csfs_LS%csfs(icsf)%subc(isubc)%iS
            S_i= csfs_LS%csfs(icsf)%iM2(isubc)
            L_i= csfs_LS%csfs(icsf)%iM1(isubc)
         end do
      end do
      if (i ==0) then
      write(iwrite_log,*)                                 &
      "There are not any configuration for evaluation"
      write(*,*)                                 &
      "There are not any configuration for evaluation"
         stop
      end if
      write(iwrite_log,*)                                 &
      "--------------------------------------------------"
      write(iwrite_log,*)                                              &
      'Total number of CSF in the input file:      ',              &
      csfs_LS%nr_of_csfs
      write(iwrite_log,*)                                              &
      'Total number of CSF included in calculation:',i
      end subroutine print_csf
!
      end module Coupling_getdata_mchf
