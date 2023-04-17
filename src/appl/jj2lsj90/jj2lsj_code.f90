!
!***********************************************************************
!                                                                      *
      MODULE jj2lsj_code
!                                                                      *
!     This module contains the procedures which are used to performe   *
!     the jj-LSJ transformation.                                       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     Vilnius                                  last update: Jan 2017   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      implicit none
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: asf2ls
!               Expands an atomic state function, which is represented
!               in a jj-coupling CSF basis into a basis of LS-coupling CSF.
      public  :: coefLSjj
!               Returns the value of the LS-jj transformation matrix
!               for a given set of quantum numbers.
      private :: coefLSjj2
!               Returns the value of the LS-jj transformation matrix
!               (l^2 LSJ| j_1 j_2 J).
      private :: coefLSjjs
!               Returns the value of the LS-jj transformation matrix
!               form an array LS_jj_*_*.
      public  :: dallocASFLS
!               Dellocates the storage of asf_set_LS.
      public  :: getchLS
!               A spectroscopic notation of shell in LS coupling is return.
      public  :: getxj
!
      public  :: gettermLS
!               This procedure return all allowed subshell terms
!               (l, w, Q, L, S) for given l^N which must be 0, 1, 2 or 3.
      public  :: inscreen
!               The input from the screen.
      public  :: inscreenlev
!               Attempts to interpret the serial level numbers from a
!               string.
      public  :: jj2lsj
!               Controls the transformation of atomic states from a jj-
!               into a LS-coupled CSF basis.
      public  :: packlsCSF
!                Encording all CSF in LS - coupling.
      public  :: prCSFjj
!               Prints all information about the single CSF in jj-coupling
      public  :: prCSFLS
!               Prints all information about the single CSF in LS-coupling
      public  :: prCSFLSall
!               Prints all information about the single CSF scheme csf_set
!               in LS-coupling
      private :: setLS
!               This subroutine fills up the variable asf_set_LS%csf_set_LS
!               with data generated using the one from asf_set%csf_set
!               .....................................................
!               This subroutine contains the following internal routines:
!                * subroutine setLS_action
!                  The subroutine defines the "action" of subroutine
!                  setLS_job_count: whether it counts
!                  the number of csfs_LS (asf_set_LS%csf_set_LS%novcsf)
!                  or fills the arrays of wave functions in LS coupling
!                  with asf_set_LS%csf_set_LS%csf(...) with
!                  the corresponding quantum nubers.
!                * subroutine setLS_add_quantum_numbers
!                  The subroutine adds quantum numbers stored
!                  in temprorary arrays Li, Si, L_i, S_i, w, Q to
!                  the corresponding arrays of asf_set_LS%csf_set_LS%csf().
!                  private :: setLS_job_count
!                * recursive subroutine setLS_job_count
!                  Recursive subroutine for the calculation of the
!                  number of csfs_LS and corresponding quantum numbers.
!                * function setLS_equivalent_csfs
!                  This subroutine defines the "equivalency" of two csfs_jj
!                  in the sence of generation of csfs_LS
!                  number of csfs_LS and corresponding quantum numbers.
!               .....................................................
      public  :: traLSjj
!               Return the value of the transformation matrix
!               from jj- to LS-coupling scheme in the case of any
!               number of open shells.
      public  :: traLSjjmp
!               Return the value of main part of the transformation
!               matrix from jj- to LS-coupling scheme in the
!               case of any number of open shells.
      private :: uniquelsj
!               Subroutine defines a unique labels for energy levels
!-----------------------------------------------
!   D e f i n i t i o n   o f   A S F
!   i n   L S -   C o u p l i n g
!-----------------------------------------------
      type, public :: nl
         integer :: n, l
      end type nl
!
      type, public :: cs_function_LS
         integer :: totalJ
         character(len=1)  :: parity
         integer, dimension(:), pointer :: occupation
         integer, dimension(:), pointer :: seniority
         integer, dimension(:), pointer :: w
         integer, dimension(:), pointer :: shellL
         integer, dimension(:), pointer :: shellS
         integer, dimension(:), pointer :: shellLX
         integer, dimension(:), pointer :: shellSX
      end type cs_function_LS
!
      type, public :: csf_basis_LS
         integer :: nocsf         ! Number of CSF in the basis.
         integer :: nwshells      ! Number of (nonrelativistic) shells.
         integer :: nwcore        ! Number of (closed) core shells.
         integer :: number_of_electrons
         type(nl), dimension(:), pointer             :: shell
         type(cs_function_LS), dimension(:), pointer :: csf
         type(parent_from_jj), dimension(:), pointer :: parent
      end type csf_basis_LS
!
      type, public :: as_function_LS
         integer          :: level_No
         integer          :: max_csf_No
         integer          :: totalL, totalS, totalJ
         character(len=1) :: parity
         real(DOUBLE)     :: energy
         real(DOUBLE), dimension(:), pointer :: eigenvector
      end type as_function_LS
!
      type, public :: asf_basis_LS
         integer      :: noasf           ! Number of considered ASF.
         real(DOUBLE) :: average_energy  ! Averaged energy of this set of ASF.
         type(as_function_LS), dimension(:), pointer :: asf
         type(csf_basis_LS) :: csf_set_LS
      end type asf_basis_LS
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      type, public :: parent_from_jj
         integer :: parent_minus
         integer :: parent_plus
      end type parent_from_jj
!
      type, public :: lsj_list
         integer::list_size                    !number items in a list
         integer, dimension(:),pointer:: items !serial numbers of lists items
      end type lsj_list
!
      type(asf_basis_LS), public :: asf_set_LS
!
      integer, private :: IMINCOMPOFF
      integer, dimension(:,:),pointer:: Jcoup ! xJ values
      real(DOUBLE) :: EPSNEW, MINCOMP
      character(len=1), dimension(0:20), private, parameter :: L_string =&
               (/ "S","P","D","F","G","H","I","K","L","M","N","O","Q",   &
                  "R","T","U","V","W","X","Y","Z" /)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: Blocks_number  = 20
      INTEGER, PARAMETER :: Vectors_number = 100000
!-----------------------------------------------
!
CONTAINS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE asf2ls(iw1,ithresh,levmax,IBLKNUM,levels,NCFMIN,     &
                        NCFMAX,ioutT)
!                                                                      *
!     Expands an atomic state functions from the same block,           *
!     which is represented in a jj-coupling CSF basis into a basis     *
!     of LS-coupling CSF.                                              *
!                                                                      *
!     Calls: ispar, itjpo, tranLSjj.                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!     Modified by G. Gaigalas                                   2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE EIGV_C,          ONLY: EVEC
      USE PRNT_C,          ONLY: IVEC
      USE ORB_C,           ONLY: NCF
      USE CONS_C,          ONLY: ZERO, ONE
!      USE JJ2LSJ_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE ispar_I
!      USE idigit_I
!      USE ittk_I
!      USE ixjtik_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: iw1, levmax,IBLKNUM,ioutT
      integer, intent(in) :: NCFMIN, NCFMAX
      integer, dimension(:), intent(in) :: ithresh
      integer, dimension(Blocks_number,Vectors_number), intent(in) :: levels
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=18) :: String
      integer      :: jj_number, lev, level, LS_number, ijj, iLS
      integer      :: LOC, IMINCOMP, IDUM
      real(DOUBLE) :: wa_transformation
      real(DOUBLE), dimension(Vectors_number) :: wa
      real(DOUBLE), dimension(Vectors_number) :: wb
!-----------------------------------------------
      wb = zero
!GGR      if(ioutT == 1) write(59) ' *   Block Number=',IBLKNUM
      if(ioutT == 1) write(59,'(A18,I6)') ' *   Block Number=',IBLKNUM
      if(ioutT == 2) then
!GGR         read(59) String, IDUM
         read(59,'(A18,I6)') String, IDUM
         if(String(1:18) /= ' *   Block Number=' .or.                  &
                                                  IDUM /= IBLKNUM) then
            print*, "Error in transformation file *.lsj.T"
         end if
      end if
      do LS_number = 1, asf_set_LS%csf_set_LS%nocsf
         if ((asf_set_LS%csf_set_LS%csf(LS_number)%parity == "+" &
            .and.  ISPAR(iw1) == 1)  .or.                        &
            (asf_set_LS%csf_set_LS%csf(LS_number)%parity  == "-" &
            .and.  ISPAR(iw1) == -1)) then
            wa = zero
            do jj_number = NCFMIN, NCFMAX
               if(ithresh(jj_number) == 1 .and.                  &
                 (asf_set_LS%csf_set_LS%csf(LS_number)%totalJ == &
                  ITJPO(jj_number)-1)) then
                  if(ioutT <= 1)                                       &
                  wa_transformation = traLSjj(jj_number,LS_number)
!GGR                  if(ioutT == 1) write(59) wa_transformation
                  if(ioutT == 1) write(59,'(2I6,2X,F16.9)')            &
                  jj_number-NCFMIN+1,LS_number,wa_transformation
                  if(ioutT == 2) then
!GGR                  read(59) wa_transformation
                     read(59,'(2I6,2X,F16.9)')                         &
                     ijj, iLS, wa_transformation
                     if(ijj /= jj_number-NCFMIN+1 .or.                 &
                        iLS /= LS_number ) then
                          print*, "Error in jj2lsj transformation file"
                          stop
                     end if
                  end if
                  do lev = 1,levmax
                     level = levels(IBLKNUM,lev)
                     LOC = (level-1)*NCF
                     wa(lev) = wa(lev)+EVEC(jj_number+LOC)*wa_transformation
                  end do
               end if
            end do
            IMINCOMP = 1
            do lev = 1,levmax
               level = levels(IBLKNUM,lev)
               wb(lev) = wa(lev)*wa(lev) + wb(lev)
               asf_set_LS%asf(level)%eigenvector(LS_number) = wa(lev)
               asf_set_LS%asf(level)%level_No = IVEC(level)
               if((wb(lev)*100.+dabs(MINCOMP)) <= 100.) IMINCOMP = 0
            end do
            if(IMINCOMPOFF == 1 ) THEN
               if(IMINCOMP == 1) GO TO 1
            end if
         end if
      end do
    1 continue
      return
      END SUBROUTINE asf2ls
!
!***********************************************************************
!                                                                      *
      FUNCTION coefLSjj(l_shell,N,w,Q,L,S,J,jm_shell,Nm,Qm,Jm,jp_shell,&
                                                     Qp,Jp) result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix for a given *
!     set of quantum numbers form an array LS_jj_*_*.                  *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers Q, L, S, J, jm_shell, Qm, Jm, jp_shell,  *
!     Qp, Jp.                                                          *
!                                                                      *
!     Calls: coefLSjj2, coefLSjjs.                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_shell, N, w, Q, L, S, J,            &
                             jm_shell, Nm, Qm, Jm, jp_shell, Qp, Jp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(DOUBLE) :: wa
      integer      :: NN, N1, N2, factor
!-----------------------------------------------
      wa = zero; factor = ONE
      if (jm_shell < jp_shell .or. jm_shell == jp_shell)  then
         if (l_shell > 0 .and. N > 2*l_shell +1)  then
            NN = 4*l_shell + 2 - N;   N1 = jm_shell + 1 - Nm
            N2 = jp_shell + 1 - N + Nm
            if (mod(2*l_shell+1-Q -((jm_shell+1)/2-Qm)         &
                 -((jp_shell+1)/2-Qp),4) /= 0) factor = - factor
         else
            NN = N;   N1 = Nm;   N2 = N - Nm;
         end if
         if (NN == 1  .or.  NN == 0)  then
            if (NN == 0 .and. N1 == 0 .and. N2 == 0) then
               wa = one
            else if (N1 == 1 .and. N2 == 0 .and. J == Jm) then
               wa = one
            else if (N1 == 0 .and. N2 == 1 .and. J == Jp) then
               wa = one
            else
               wa = zero
            end if
         else if (NN == 2)  then
            if (J > Jm + Jp  .or.  J < abs(Jm - Jp)) then
               wa = zero
            else
               if (N1 == 2 .and. N2 == 0) then
                  wa = coefLSjj2(l_shell,L,S,J,jm_shell,jm_shell)
               else if (N1 == 1 .and. N2 == 1) then
                  wa = coefLSjj2(l_shell,L,S,J,jm_shell,jp_shell)
               else if (N1 == 0 .and. N2 == 2) then
                  wa = coefLSjj2(l_shell,L,S,J,jp_shell,jp_shell)
               end if;
            end if
         else if (l_shell==1 .or. l_shell==2 .or. l_shell==3)  then
            wa = coefLSjjs(l_shell,NN,w,Q,L,S,J,N1,Qm,Jm,Qp,Jp)
         end if
      else if (l == 0)  then
         if (S == J .and. jm_shell == 1 .and. NN == N1)  then
            wa = one
         else if (S == J .and. jp_shell == 1 .and. NN == N2)  then
            wa = one
         else
            wa = zero
         end if
      else
         stop "coefLSjj: program stop A."
      end if
      wa = factor * wa
      END FUNCTION coefLSjj
!
!***********************************************************************
!                                                                      *
      FUNCTION coefLSjj2(l_shell,L,S,J,jm_shell,jp_shell)   result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix             *
!                                                (l^2 LSJ| j_1 j_2 J). *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers L, S, J, jm_shell, jp_shell.             *
!                                                                      *
!     Calls: nine.                                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_shell, L, S, J, jm_shell, jp_shell
      real(DOUBLE)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: delta_J
      real(DOUBLE) :: RAC9
!-----------------------------------------------
      wa = zero
      if (mod(L+S,4) /= 0)  then
         wa = zero
      elseif (mod(L+S+J,2) /= 0)  then
         wa = zero
      elseif (mod(jm_shell+jp_shell+J,2) /= 0)  then
         wa = zero
      else
         if (jm_shell == jp_shell) then
            if (mod(J,4) /= 0)  then
               wa = zero
            else
               wa = (jm_shell+one) * sqrt((L+one) * (S+one))
            end if
         else
            wa = sqrt(2*(L+one)*(S+one)*(jm_shell+one)*(jp_shell+one))
         end if
         call nine(2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J,   &
                   1,delta_J,RAC9)
         if (delta_J /= 0)  then
            call nine(2*l_shell,2*l_shell,L,1,1,S,jm_shell,jp_shell,J,&
                      0,delta_J,RAC9)
            wa = wa * RAC9
         end if
      end if
      END FUNCTION coefLSjj2
!
!***********************************************************************
!                                                                      *
      FUNCTION coefLSjjs(lshell,N,w,Q,L,S,J,Nm,Qm,Jm,Qp,Jp)  result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix for a given *
!     set of quantum numbers.                                          *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers Q, L, S, J, Qm, Jm, Qp, Jp.              *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE jj2lsj_C
      USE jj2lsj_data_1_C
      USE jj2lsj_data_2_C
      USE jj2lsj_data_3_C
      USE CONS_C,          ONLY: ZERO, ONE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: lshell, N, w, Q, L, S, J, Nm, Qm, Jm, Qp, Jp
      real(DOUBLE)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      wa = zero
      if (lshell == 0) then
      else if (lshell == 1) then
         select case(N)
         case(3)
!     Use data from the array LS_jj_p_3
            do  i = 1,LS_jj_number_p3
               if(w  ==LS_jj_p_3(i)%w  .and. Q  ==LS_jj_p_3(i)%Q  .and. &
                  L  ==LS_jj_p_3(i)%L  .and. S  ==LS_jj_p_3(i)%S  .and. &
                  J  ==LS_jj_p_3(i)%J  .and. Nm ==LS_jj_p_3(i)%Nm .and. &
                  Qm ==LS_jj_p_3(i)%Qm .and. Jm ==LS_jj_p_3(i)%Jm .and. &
                  Qp ==LS_jj_p_3(i)%Qp .and. Jp ==LS_jj_p_3(i)%Jp) then
                  wa = LS_jj_p_3(i)%factor * &
                       sqrt( one*LS_jj_p_3(i)%nom/LS_jj_p_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_p_4
            do  i = 1,LS_jj_number_p4
               if(w  ==LS_jj_p_4(i)%w  .and. Q  ==LS_jj_p_4(i)%Q  .and. &
                  L  ==LS_jj_p_4(i)%L  .and. S  ==LS_jj_p_4(i)%S  .and. &
                  J  ==LS_jj_p_4(i)%J  .and. Nm ==LS_jj_p_4(i)%Nm .and. &
                  Qm ==LS_jj_p_4(i)%Qm .and. Jm ==LS_jj_p_4(i)%Jm .and. &
                  Qp ==LS_jj_p_4(i)%Qp .and. Jp ==LS_jj_p_4(i)%Jp) then
                  wa = LS_jj_p_4(i)%factor * &
                       sqrt( one*LS_jj_p_4(i)%nom/LS_jj_p_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_p_5
            do  i = 1,LS_jj_number_p5
               if(w  ==LS_jj_p_5(i)%w  .and. Q  ==LS_jj_p_5(i)%Q  .and. &
                  L  ==LS_jj_p_5(i)%L  .and. S  ==LS_jj_p_5(i)%S  .and. &
                  J  ==LS_jj_p_5(i)%J  .and. Nm ==LS_jj_p_5(i)%Nm .and. &
                  Qm ==LS_jj_p_5(i)%Qm .and. Jm ==LS_jj_p_5(i)%Jm .and. &
                  Qp ==LS_jj_p_5(i)%Qp .and. Jp ==LS_jj_p_5(i)%Jp) then
                  wa = LS_jj_p_5(i)%factor * &
                       sqrt( one*LS_jj_p_5(i)%nom/LS_jj_p_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_p_6
            do  i = 1,LS_jj_number_p6
               if(w  ==LS_jj_p_6(i)%w  .and. Q  ==LS_jj_p_6(i)%Q  .and. &
                  L  ==LS_jj_p_6(i)%L  .and. S  ==LS_jj_p_6(i)%S  .and. &
                  J  ==LS_jj_p_6(i)%J  .and. Nm ==LS_jj_p_6(i)%Nm .and. &
                  Qm ==LS_jj_p_6(i)%Qm .and. Jm ==LS_jj_p_6(i)%Jm .and. &
                  Qp ==LS_jj_p_6(i)%Qp .and. Jp ==LS_jj_p_6(i)%Jp) then
                  wa = LS_jj_p_6(i)%factor * &
                       sqrt( one*LS_jj_p_6(i)%nom/LS_jj_p_6(i)%denom )
                  return
               end if
            end do
         case default
            stop "coefLSjjs: program stop A."
         end select
      else if (lshell == 2) then
         select case(N)
         case(3)
!     Use data from the array LS_jj_d_3
            do  i = 1,LS_jj_number_d3
               if(w  ==LS_jj_d_3(i)%w  .and. Q  ==LS_jj_d_3(i)%Q  .and. &
                  L  ==LS_jj_d_3(i)%L  .and. S  ==LS_jj_d_3(i)%S  .and. &
                  J  ==LS_jj_d_3(i)%J  .and. Nm ==LS_jj_d_3(i)%Nm .and. &
                  Qm ==LS_jj_d_3(i)%Qm .and. Jm ==LS_jj_d_3(i)%Jm .and. &
                  Qp ==LS_jj_d_3(i)%Qp .and. Jp ==LS_jj_d_3(i)%Jp) then
                  wa = LS_jj_d_3(i)%factor * &
                       sqrt( one*LS_jj_d_3(i)%nom/LS_jj_d_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_d_4
            do  i = 1,LS_jj_number_d4
               if(w  ==LS_jj_d_4(i)%w  .and. Q  ==LS_jj_d_4(i)%Q  .and. &
                  L  ==LS_jj_d_4(i)%L  .and. S  ==LS_jj_d_4(i)%S  .and. &
                  J  ==LS_jj_d_4(i)%J  .and. Nm ==LS_jj_d_4(i)%Nm .and. &
                  Qm ==LS_jj_d_4(i)%Qm .and. Jm ==LS_jj_d_4(i)%Jm .and. &
                  Qp ==LS_jj_d_4(i)%Qp .and. Jp ==LS_jj_d_4(i)%Jp) then
                  wa = LS_jj_d_4(i)%factor * &
                       sqrt( one*LS_jj_d_4(i)%nom/LS_jj_d_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_d_5
            do  i = 1,LS_jj_number_d5
               if(w  ==LS_jj_d_5(i)%w  .and. Q  ==LS_jj_d_5(i)%Q  .and. &
                  L  ==LS_jj_d_5(i)%L  .and. S  ==LS_jj_d_5(i)%S  .and. &
                  J  ==LS_jj_d_5(i)%J  .and. Nm ==LS_jj_d_5(i)%Nm .and. &
                  Qm ==LS_jj_d_5(i)%Qm .and. Jm ==LS_jj_d_5(i)%Jm .and. &
                  Qp ==LS_jj_d_5(i)%Qp .and. Jp ==LS_jj_d_5(i)%Jp) then
                  wa = LS_jj_d_5(i)%factor * &
                       sqrt( one*LS_jj_d_5(i)%nom/LS_jj_d_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_d_6
            do  i = 1,LS_jj_number_d6
               if(w  ==LS_jj_d_6(i)%w  .and. Q  ==LS_jj_d_6(i)%Q  .and. &
                  L  ==LS_jj_d_6(i)%L  .and. S  ==LS_jj_d_6(i)%S  .and. &
                  J  ==LS_jj_d_6(i)%J  .and. Nm ==LS_jj_d_6(i)%Nm .and. &
                  Qm ==LS_jj_d_6(i)%Qm .and. Jm ==LS_jj_d_6(i)%Jm .and. &
                  Qp ==LS_jj_d_6(i)%Qp .and. Jp ==LS_jj_d_6(i)%Jp) then
                  wa = LS_jj_d_6(i)%factor * &
                       sqrt( one*LS_jj_d_6(i)%nom/LS_jj_d_6(i)%denom )
                  return
               end if
            end do
         case(7)
!     Use data from the array LS_jj_d_7
            do  i = 1,LS_jj_number_d7
               if(w  ==LS_jj_d_7(i)%w  .and. Q  ==LS_jj_d_7(i)%Q  .and. &
                  L  ==LS_jj_d_7(i)%L  .and. S  ==LS_jj_d_7(i)%S  .and. &
                  J  ==LS_jj_d_7(i)%J  .and. Nm ==LS_jj_d_7(i)%Nm .and. &
                  Qm ==LS_jj_d_7(i)%Qm .and. Jm ==LS_jj_d_7(i)%Jm .and. &
                  Qp ==LS_jj_d_7(i)%Qp .and. Jp ==LS_jj_d_7(i)%Jp) then
                  wa = LS_jj_d_7(i)%factor * &
                       sqrt( one*LS_jj_d_7(i)%nom/LS_jj_d_7(i)%denom )
                  return
               end if
            end do
         case(8)
!     Use data from the array LS_jj_d_8
            do  i = 1,LS_jj_number_d8
               if(w  ==LS_jj_d_8(i)%w  .and. Q  ==LS_jj_d_8(i)%Q  .and. &
                  L  ==LS_jj_d_8(i)%L  .and. S  ==LS_jj_d_8(i)%S  .and. &
                  J  ==LS_jj_d_8(i)%J  .and. Nm ==LS_jj_d_8(i)%Nm .and. &
                  Qm ==LS_jj_d_8(i)%Qm .and. Jm ==LS_jj_d_8(i)%Jm .and. &
                  Qp ==LS_jj_d_8(i)%Qp .and. Jp ==LS_jj_d_8(i)%Jp) then
                  wa = LS_jj_d_8(i)%factor * &
                       sqrt( one*LS_jj_d_8(i)%nom/LS_jj_d_8(i)%denom )
                  return
               end if
            end do
         case(9)
!     Use data from the array LS_jj_d_9
            do  i = 1,LS_jj_number_d9
               if(w  ==LS_jj_d_9(i)%w  .and. Q  ==LS_jj_d_9(i)%Q  .and. &
                  L  ==LS_jj_d_9(i)%L  .and. S  ==LS_jj_d_9(i)%S  .and. &
                  J  ==LS_jj_d_9(i)%J  .and. Nm ==LS_jj_d_9(i)%Nm .and. &
                  Qm ==LS_jj_d_9(i)%Qm .and. Jm ==LS_jj_d_9(i)%Jm .and. &
                  Qp ==LS_jj_d_9(i)%Qp .and. Jp ==LS_jj_d_9(i)%Jp) then
                  wa = LS_jj_d_9(i)%factor * &
                       sqrt( one*LS_jj_d_9(i)%nom/LS_jj_d_9(i)%denom )
                  return
               end if
            end do
         case(10)
!     Use data from the array LS_jj_d_10
            do  i = 1,LS_jj_number_d10
               if(w ==LS_jj_d_10(i)%w  .and. Q ==LS_jj_d_10(i)%Q  .and. &
                  L ==LS_jj_d_10(i)%L  .and. S ==LS_jj_d_10(i)%S  .and. &
                  J ==LS_jj_d_10(i)%J  .and. Nm==LS_jj_d_10(i)%Nm .and. &
                  Qm==LS_jj_d_10(i)%Qm .and. Jm==LS_jj_d_10(i)%Jm .and. &
                  Qp==LS_jj_d_10(i)%Qp .and. Jp==LS_jj_d_10(i)%Jp) then
                  wa= LS_jj_d_10(i)%factor * &
                       sqrt( one*LS_jj_d_10(i)%nom/LS_jj_d_10(i)%denom )
                  return
               end if
            end do
         case default
            stop "coefLSjjs: program stop B."
         end select
      else if (lshell == 3) then
         select case(N)
         case(3)
!     Use data from the array LS_jj_f_3
            do  i = 1,LS_jj_number_f3
               if(w  ==LS_jj_f_3(i)%w  .and. Q  ==LS_jj_f_3(i)%Q  .and. &
                  L  ==LS_jj_f_3(i)%L  .and. S  ==LS_jj_f_3(i)%S  .and. &
                  J  ==LS_jj_f_3(i)%J  .and. Nm ==LS_jj_f_3(i)%Nm .and. &
                  Qm ==LS_jj_f_3(i)%Qm .and. Jm ==LS_jj_f_3(i)%Jm .and. &
                  Qp ==LS_jj_f_3(i)%Qp .and. Jp ==LS_jj_f_3(i)%Jp) then
                  wa = LS_jj_f_3(i)%factor * &
                       sqrt( one*LS_jj_f_3(i)%nom/LS_jj_f_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_f_4
            do  i = 1,LS_jj_number_f4
               if(w  ==LS_jj_f_4(i)%w  .and. Q  ==LS_jj_f_4(i)%Q  .and. &
                  L  ==LS_jj_f_4(i)%L  .and. S  ==LS_jj_f_4(i)%S  .and. &
                  J  ==LS_jj_f_4(i)%J  .and. Nm ==LS_jj_f_4(i)%Nm .and. &
                  Qm ==LS_jj_f_4(i)%Qm .and. Jm ==LS_jj_f_4(i)%Jm .and. &
                  Qp ==LS_jj_f_4(i)%Qp .and. Jp ==LS_jj_f_4(i)%Jp) then
                  wa = LS_jj_f_4(i)%factor * &
                       sqrt( one*LS_jj_f_4(i)%nom/LS_jj_f_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_f_5
            do  i = 1,LS_jj_number_f5
               if(w  ==LS_jj_f_5(i)%w  .and. Q  ==LS_jj_f_5(i)%Q  .and. &
                  L  ==LS_jj_f_5(i)%L  .and. S  ==LS_jj_f_5(i)%S  .and. &
                  J  ==LS_jj_f_5(i)%J  .and. Nm ==LS_jj_f_5(i)%Nm .and. &
                  Qm ==LS_jj_f_5(i)%Qm .and. Jm ==LS_jj_f_5(i)%Jm .and. &
                  Qp ==LS_jj_f_5(i)%Qp .and. Jp ==LS_jj_f_5(i)%Jp) then
                  wa = LS_jj_f_5(i)%factor * &
                       sqrt( one*LS_jj_f_5(i)%nom/LS_jj_f_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_f_6
            do  i = 1,LS_jj_number_f6
               if(w  ==LS_jj_f_6(i)%w  .and. Q  ==LS_jj_f_6(i)%Q  .and. &
                  L  ==LS_jj_f_6(i)%L  .and. S  ==LS_jj_f_6(i)%S  .and. &
                  J  ==LS_jj_f_6(i)%J  .and. Nm ==LS_jj_f_6(i)%Nm .and. &
                  Qm ==LS_jj_f_6(i)%Qm .and. Jm ==LS_jj_f_6(i)%Jm .and. &
                  Qp ==LS_jj_f_6(i)%Qp .and. Jp ==LS_jj_f_6(i)%Jp) then
                  wa = LS_jj_f_6(i)%factor * &
                       sqrt( one*LS_jj_f_6(i)%nom/LS_jj_f_6(i)%denom )
                  return
               end if
            end do
         case(7)
!     Use data from the array LS_jj_f_7
            do  i = 1, LS_jj_number_f7, 1
               if (w == LS_jj_f_7(i)%w) then
                  if(Q == LS_jj_f_7(i)%Q) then
                     if(L == LS_jj_f_7(i)%L) then
                        if(S == LS_jj_f_7(i)%S) then
                           if(J == LS_jj_f_7(i)%J) then
                              if(Nm == LS_jj_f_7(i)%Nm) then
                                 if(Qm == LS_jj_f_7(i)%Qm .and. &
                                    Jm == LS_jj_f_7(i)%Jm .and. &
                                    Qp == LS_jj_f_7(i)%Qp .and. &
                                    Jp == LS_jj_f_7(i)%Jp) then
                                       wa = one * LS_jj_f_7(i)%factor  * &
                                            dsqrt(one*LS_jj_f_7(i)%nom)/ &
                                            dsqrt(one*LS_jj_f_7(i)%denom )
                                       return
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
            end do
         case default
            stop "coefLSjjs: program stop C."
         end select
      else
         stop "coefLSjjs: program stop D."
      end if
      END FUNCTION coefLSjjs
!
!***********************************************************************
!                                                                      *
      SUBROUTINE dallocASFLS(asf_set_LS)
!                                                                      *
!     Dellocates the storage of asf_set_LS.                            *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,    ONLY: DOUBLE
      USE PRNT_C,             ONLY: NVEC
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      type(asf_basis_LS), intent(inout) :: asf_set_LS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      deallocate(asf_set_LS%csf_set_LS%shell)
      deallocate(asf_set_LS%csf_set_LS%parent)
      do i = 1, asf_set_LS%csf_set_LS%nocsf, 1
         deallocate(asf_set_LS%csf_set_LS%csf(i)%occupation)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%seniority)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellL)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellS)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellLX)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellSX)
      end do
      deallocate(asf_set_LS%csf_set_LS%csf)
      do i = 1, NVEC
         deallocate(asf_set_LS%asf(i)%eigenvector)
      end do
      deallocate(asf_set_LS%asf)
      END SUBROUTINE dallocASFLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE getchLS(csf_number,shell_number,string_LS, string_XLS)
!                                                                      *
!     The spectroscopic notation of a shell in LS coupling is returned.*
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: csf_number,shell_number
      character(len=4), intent(out) :: string_LS, string_XLS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: string_lenth
      character(len=1) :: string_S, string_v
      character(len=4) :: string_CNUM
!-----------------------------------------------
      call                                                            &
      convrt_double(2*(1+asf_set_LS%csf_set_LS%csf(csf_number)%shellS(shell_number)),&
      string_CNUM,string_lenth)
      string_S = string_CNUM(1:string_lenth)
      if (asf_set_LS%csf_set_LS%shell(shell_number)%l < 3)    then
         call                                                         &
         convrt_double(2*asf_set_LS%csf_set_LS%csf(csf_number)%seniority(shell_number),&
         string_CNUM,string_lenth)
         string_v = string_CNUM(1:string_lenth)
      else
         call                                                         &
         convrt_double(2*asf_set_LS%csf_set_LS%csf(csf_number)%w(shell_number),&
         string_CNUM,string_lenth)
         string_v = string_CNUM(1:string_lenth)
      end if
      string_LS = string_S //                                         &
      L_string(asf_set_LS%csf_set_LS%csf(csf_number)%shellL(shell_number)/2)&
      //string_v
      !
      call                                                            &
      convrt_double(2*(1+asf_set_LS%csf_set_LS%csf(csf_number)%shellSX(shell_number)),&
      string_CNUM,string_lenth)
      string_S =  string_CNUM(1:string_lenth)
      string_XLS = string_S //                                        &
         L_string(asf_set_LS%csf_set_LS%csf(csf_number)%shellLX(shell_number)/2)
      END SUBROUTINE getchLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE getxj
!                                                                      *
!                                                                      *
!     Calls: jcup, jqs, ichop.                                         *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE M_C,             ONLY: NCORE
      USE ORB_C,           ONLY: NCF, NW
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE jqs_I
      USE ichop_I
      USE jcup_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JCNT, JCNTOP, JNCF, JNW
!-----------------------------------------------
      Jcoup = 0
      DO JNCF = 1, NCF
         JCNT = 1
         JCNTOP = 0
         DO JNW = NCORE+1, NW
            IF(JNW == 1) THEN
               IF(ICHOP(JNW,JNCF) /= 0) THEN
                  IF(NW == 1) THEN
                      Jcoup(JNW,JNCF) = 0
                  ELSE
                      IF(ICHOP(JNW+1,JNCF) /= 0) THEN
                         Jcoup(JNW,JNCF)   = 0
                         Jcoup(JNW+1,JNCF) = 0
                      ELSE
                         Jcoup(JNW,JNCF)   = 0
                         Jcoup(JNW+1,JNCF) = JQS(3,JNW+1,JNCF) - 1
                         JCNTOP = 1
                     END IF
                  END IF
               ELSE
                  JCNTOP = 1
                  IF(NW > 1) THEN
                     IF (ICHOP(JNW+1,JNCF) == 0) THEN
!GG 2015_06_21 Gediminas Gaigalas
!GG                        Jcoup(JNW,JNCF)   = JCUP(JCNT,JNCF) - 1
                        Jcoup(JNW,JNCF)   = JQS(3,JNW,JNCF) - 1
                        Jcoup(JNW+1,JNCF) = JCUP(JCNT,JNCF) - 1
                        JCNT = JCNT + 1
                     ELSE
                        Jcoup(JNW,JNCF)   = JQS(3,JNW,JNCF) - 1
                        Jcoup(JNW+1,JNCF) = JQS(3,JNW,JNCF) - 1
                     ENDIF
                  ELSE
                     Jcoup(JNW,JNCF) = JCUP(JCNT,JNCF) - 1
                  ENDIF
               ENDIF
            ELSE IF(JNW == 2 .AND. NCORE+1 .EQ. 2) THEN
               IF(ICHOP(JNW,JNCF) /= 0) THEN
                  Jcoup(JNW,JNCF) = 0
               ELSE
                  JCNTOP = 1
                  Jcoup(JNW,JNCF)   = JQS(3,JNW,JNCF) - 1
               ENDIF
            ELSE IF(JNW > 2) THEN
               IF (ICHOP(JNW,JNCF) /= 0) THEN
                  IF(JNW == NCORE+1) THEN
                     Jcoup(JNW,JNCF) = JQS(3,JNW,JNCF) - 1
                  ELSE
                     Jcoup(JNW,JNCF) = Jcoup(JNW-1,JNCF)
                  END IF
               ELSE
                  IF (JCNTOP /= 0) THEN
                     Jcoup(JNW,JNCF) = JCUP(JCNT,JNCF)  - 1
                     JCNT = JCNT + 1
                  ELSE
                     Jcoup(JNW,JNCF) = JQS(3,JNW,JNCF) - 1
                  ENDIF
                  JCNTOP = JCNTOP + 1
               END IF
            END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE  getxj
!
!***********************************************************************
!                                                                      *
      SUBROUTINE gettermLS(l_shell,N,LS,number)
!                                                                      *
!     This procedure returns all allowed subshell terms (l,w,Q,L,S)    *
!     for given l^N, N = 0, 1, 2 or 3.                                 *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE jj2lsj_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)                                 :: l_shell, N
      type(subshell_term_LS), dimension(120), intent(out) :: LS
      integer, intent(out)                                :: number
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: M_Q, i, j
!-----------------------------------------------
      M_Q = N - 2* l_shell - 1;  j = 0
      select case (l_shell)
      case (0)
         do i = 1,2
            if (mod(M_Q + term_LS_s(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_s(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_s(i)%l_shell
                  LS(j)%w        = term_LS_s(i)%w
                  LS(j)%Q        = term_LS_s(i)%Q
                  LS(j)%LL       = term_LS_s(i)%LL
                  LS(j)%S        = term_LS_s(i)%S
               end if
            end if
         end do
      case (1)
         do i = 1,6
            if (mod(M_Q + term_LS_p(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_p(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_p(i)%l_shell
                  LS(j)%w        = term_LS_p(i)%w
                  LS(j)%Q        = term_LS_p(i)%Q
                  LS(j)%LL       = term_LS_p(i)%LL
                  LS(j)%S        = term_LS_p(i)%S
               end if
            end if
         end do
      case (2)
         do i = 1,32
            if (mod(M_Q + term_LS_d(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_d(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_d(i)%l_shell
                  LS(j)%w        = term_LS_d(i)%w
                  LS(j)%Q        = term_LS_d(i)%Q
                  LS(j)%LL       = term_LS_d(i)%LL
                  LS(j)%S        = term_LS_d(i)%S
               end if
            end if
         end do
      case (3)
         do i = 1,238
            if (mod(M_Q + term_LS_f(i)%Q,2) == 0) then
               if (abs(M_Q) <= term_LS_f(i)%Q) then
                  j = j + 1
                  LS(j)%l_shell  = term_LS_f(i)%l_shell
                  LS(j)%w        = term_LS_f(i)%w
                  LS(j)%Q        = term_LS_f(i)%Q
                  LS(j)%LL       = term_LS_f(i)%LL
                  LS(j)%S        = term_LS_f(i)%S
               end if
            end if
         end do
      case (4)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_g1(i)%l_shell
            LS(j)%w        = term_LS_g1(i)%w
            LS(j)%Q        = term_LS_g1(i)%Q
            LS(j)%LL       = term_LS_g1(i)%LL
            LS(j)%S        = term_LS_g1(i)%S
         case (2)
            do i = 1,9
               j = j + 1
               LS(j)%l_shell  = term_LS_g2(i)%l_shell
               LS(j)%w        = term_LS_g2(i)%w
               LS(j)%Q        = term_LS_g2(i)%Q
               LS(j)%LL       = term_LS_g2(i)%LL
               LS(j)%S        = term_LS_g2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop A."
         end select
      case (5)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_h1(i)%l_shell
            LS(j)%w        = term_LS_h1(i)%w
            LS(j)%Q        = term_LS_h1(i)%Q
            LS(j)%LL       = term_LS_h1(i)%LL
            LS(j)%S        = term_LS_h1(i)%S
         case (2)
            do i = 1,11
               j = j + 1
               LS(j)%l_shell  = term_LS_h2(i)%l_shell
               LS(j)%w        = term_LS_h2(i)%w
               LS(j)%Q        = term_LS_h2(i)%Q
               LS(j)%LL       = term_LS_h2(i)%LL
               LS(j)%S        = term_LS_h2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop B."
         end select
      case (6)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_i1(i)%l_shell
            LS(j)%w        = term_LS_i1(i)%w
            LS(j)%Q        = term_LS_i1(i)%Q
            LS(j)%LL       = term_LS_i1(i)%LL
            LS(j)%S        = term_LS_i1(i)%S
         case (2)
            do i = 1,13
               j = j + 1
               LS(j)%l_shell  = term_LS_i2(i)%l_shell
               LS(j)%w        = term_LS_i2(i)%w
               LS(j)%Q        = term_LS_i2(i)%Q
               LS(j)%LL       = term_LS_i2(i)%LL
               LS(j)%S        = term_LS_i2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop C."
         end select
      case (7)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_k1(i)%l_shell
            LS(j)%w        = term_LS_k1(i)%w
            LS(j)%Q        = term_LS_k1(i)%Q
            LS(j)%LL       = term_LS_k1(i)%LL
            LS(j)%S        = term_LS_k1(i)%S
         case (2)
            do i = 1,15
               j = j + 1
               LS(j)%l_shell  = term_LS_k2(i)%l_shell
               LS(j)%w        = term_LS_k2(i)%w
               LS(j)%Q        = term_LS_k2(i)%Q
               LS(j)%LL       = term_LS_k2(i)%LL
               LS(j)%S        = term_LS_k2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop D."
         end select
      case (8)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_l1(i)%l_shell
            LS(j)%w        = term_LS_l1(i)%w
            LS(j)%Q        = term_LS_l1(i)%Q
            LS(j)%LL       = term_LS_l1(i)%LL
            LS(j)%S        = term_LS_l1(i)%S
         case (2)
            do i = 1,17
               j = j + 1
               LS(j)%l_shell  = term_LS_l2(i)%l_shell
               LS(j)%w        = term_LS_l2(i)%w
               LS(j)%Q        = term_LS_l2(i)%Q
               LS(j)%LL       = term_LS_l2(i)%LL
               LS(j)%S        = term_LS_l2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop E."
         end select
      case (9)
         select case (N)
         case (1)
            i = 1;       j = 1
            LS(j)%l_shell  = term_LS_m1(i)%l_shell
            LS(j)%w        = term_LS_m1(i)%w
            LS(j)%Q        = term_LS_m1(i)%Q
            LS(j)%LL       = term_LS_m1(i)%LL
            LS(j)%S        = term_LS_m1(i)%S
         case (2)
            do i = 1,19
               j = j + 1
               LS(j)%l_shell  = term_LS_m2(i)%l_shell
               LS(j)%w        = term_LS_m2(i)%w
               LS(j)%Q        = term_LS_m2(i)%Q
               LS(j)%LL       = term_LS_m2(i)%LL
               LS(j)%S        = term_LS_m2(i)%S
            end do
         case default
            stop  "gettermLS(): program stop F."
         end select
      case default
         stop  "gettermLS(): program stop G."
      end select
      number = j
      END SUBROUTINE gettermLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE inscreen(THRESH,levels,number_of_levels,ioutC,ioutj, &
                                                          UNIQUE,ioutT)
!                                                                      *
!     The input from the screen.                                       *
!                                                                      *
!     Calls: sercsla, getxj, getmixblock, inscreenlev, openfl.         *
!     inscreenlev.                                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!     Modified by G. Gaigalas                                   2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE ORB_C,           ONLY: NCF, NW
      USE PRNT_C,          ONLY: NVEC
      USE IOUNIT_C,        ONLY: ISTDI, ISTDE
      USE CONS_C,          ONLY: EPS, ZERO
!GG      USE BLK_C,           ONLY: NEVINBLK, NBLOCK
      USE BLK_C,           ONLY: NEVINBLK, NCFINBLK, NBLOCK, TWO_J
      USE m_C,             ONLY: NCORE
      USE def_C,           ONLY: Z, NELEC
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out)      :: ioutC,ioutj,UNIQUE,ioutT
      real(DOUBLE), intent(out) :: THRESH
      integer,  dimension(Blocks_number), intent(out) :: number_of_levels
      integer, dimension(Blocks_number,Vectors_number), intent(out) :: levels
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer            :: NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
      integer            :: JBDUM, NCFINBLKDUM, NEVINBLKDUM, IATJPDUM
      integer            :: I, II, ISUM, K, NCI, ierr, JB
      integer            :: IBlock, number_of_levels_tmp
      logical            :: yes, fail, GETYN
      CHARACTER(LEN=6)   :: G92MIX
      character(len=24)  :: NAME
      character(len=256) :: record, util_csl_file
      integer, dimension(Blocks_number)    :: posi
      integer, dimension(1:Vectors_number) :: levels_tmp
!-----------------------------------------------
      NBLOCK = 0
    1 WRITE (ISTDE,*) 'Name of state'
      READ(ISTDI,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE (istde,*) 'Names may not start with a blank'
         GOTO 1
      ENDIF
!
!   Open, check, load data from, and close, the  .csl  file
      CALL SETCSLA(NAME,NCORE)
      allocate(Jcoup(1:NW,1:NCF))
      CALL GETXJ
      WRITE (ISTDE,*)
      WRITE (ISTDE,*) 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
!GG-2017
      WRITE (ISTDE,*) 'Do you need a unique labeling? (y/n)'
      YES = GETYN ()
      IF (YES) THEN
         UNIQUE = 1
      ELSE
         UNIQUE = 0
      ENDIF
!GG-2017 end
!
!   Get the eigenvectors
      CALL GETMIXBLOCK (NAME, NCI)
      WRITE (ISTDE,*) 'Default settings?  (y/n) '
      YES = GETYN ()
      IF (YES) THEN
         EPSNEW = 0.005D00
         ISUM = 0
         IMINCOMPOFF = 1
         MINCOMP = 1
         THRESH = 0.001D00
         ioutC = 0
         ioutj = 0
         ioutT = 0
         DO I = 1, NBLOCK
            number_of_levels(I) = NEVINBLK(I)
            IF(NEVINBLK(I) /= 0) THEN
               DO II = 1, NEVINBLK(I)
                  ISUM = ISUM + 1
                  levels(I,II) = ISUM
               END DO
           END IF
         END DO
      ELSE
         WRITE (ISTDE,*) 'All levels (Y/N)'
         YES = GETYN ()
         WRITE (istde,*)
         IF (YES) THEN
            ISUM = 0
            DO I = 1, NBLOCK
                number_of_levels(I) = NEVINBLK(I)
                IF(NEVINBLK(I) /= 0) THEN
                   DO II = 1, NEVINBLK(I)
                      ISUM = ISUM + 1
                      levels(I,II) = ISUM
                   END DO
               END IF
            END DO
         ELSE
            number_of_levels = 0
            posi(1) = 0
            DO I = 1, NBLOCK-1
               posi(I+1) = NEVINBLK(I) + posi(I)
            END DO
            WRITE (ISTDE,*) "Maximum number of ASFs is:",NVEC
            WRITE (ISTDE,*) "Enter the level numbers of the ASF which are to be transformed,"
    2       WRITE (ISTDE,*) "Enter the block number"
            read (*, "(I2)") IBLOCK
            WRITE (ISTDE,*) "The block number is:",IBLOCK
            WRITE (ISTDE,*) " e.g. 1 3 4  7 - 20  48  69 - 85 :"
            read (*, "(a)") record
            call inscreenlev(record,levels_tmp,number_of_levels_tmp,fail)
            if (fail) then
               WRITE (ISTDE,*) "Unable to interprete the serial level numbers; redo ..."
               goto 2
            end if
            number_of_levels(IBLOCK) = number_of_levels_tmp
            if (NVEC < number_of_levels(IBLOCK)) then
               WRITE (ISTDE,*) "There are to much ASF:", number_of_levels(IBLOCK)
               go to 2
            end if
            DO I = 1,number_of_levels_tmp
               levels(IBLOCK,I) = posi(IBLOCK) + levels_tmp(I)
            END DO
            WRITE (ISTDE,*) ""
            WRITE (ISTDE,*) "Do you need to include more levels?  (y/n)"
            YES = GETYN ()
            IF (YES) go to 2
         END IF
    3    WRITE (ISTDE,*) 'Maximum % of omitted composition'
         READ *, MINCOMP
         IF(MINCOMP == ZERO) THEN
            IMINCOMPOFF = 0
            EPSNEW = EPS*EPS
            ioutC = 1
            ioutj = 1
         ELSE IF( MINCOMP > ZERO) THEN
            IMINCOMPOFF = 1
            WRITE (ISTDE,*) 'What is the value below which an eigenvector component'
            WRITE (ISTDE,*) 'is to be neglected in the determination of the LSJ expansion:'
            WRITE (ISTDE,'(A,F8.5)') ' should be smaller than:',MINCOMP*0.01
            READ *, EPSNEW
         ELSE
            WRITE (ISTDE,*) " The maximum of omitted composition can be 100%"
            GO TO 3
         END IF
         WRITE (ISTDE,*) 'What is the value below which an eigenvector composition'
         WRITE (ISTDE,*) 'is to be neglected for printing?'
         READ *, THRESH
         IF( MINCOMP > ZERO) THEN
            WRITE (ISTDE,*) "Do you need the output file *.lsj.c?  (y/n)"
            YES = GETYN ()
            IF (YES) THEN
               ioutC = 1
            ELSE
               ioutC = 0
            END IF
            WRITE (ISTDE,*) "Do you need the output file *.lsj.j?  (y/n)"
            YES = GETYN ()
            IF (YES) THEN
               ioutj = 1
            ELSE
               ioutj = 0
            END IF
         END IF
         WRITE (ISTDE,*)                                              &
         "Do you need the transformation output file *.lsj.T?  (y/n)"
         YES = GETYN ()
         IF (YES) THEN
            ioutT = 1
         ELSE
            WRITE (ISTDE,*)                                           &
            "Will you use the transformation file *.lsj.T?  (y/n)"
            YES = GETYN ()
            IF (YES) THEN
               ioutT = 2
            ELSE
               ioutT = 0
            END IF
         END IF
      ENDIF
!
      WRITE (ISTDE,*)
      IF (IMINCOMPOFF == 1)             &
         WRITE (*,"(A,F8.3)") ' Maximum % of omitted composition is ',MINCOMP
      WRITE (*,"(A,ES8.1,A)")   &
        ' Below ',EPSNEW,' the eigenvector component is to be neglected for calculating'
      WRITE (*,"(A,ES8.1,A)")   &
        ' Below ',THRESH,' the eigenvector composition is to be neglected for printing'
      print *, " "
!
!     Opening the files  *.lsj.c
!
      IF(ioutC == 1) THEN
         util_csl_file = NAME(1:K-1)//'.lsj'//'.c'
         call OPENFL(56,util_csl_file,"formatted","UNKNOWN",ierr)
         IF (IERR .EQ. 1) THEN
            print *, 'Error when opening',util_csl_file
            STOP
         ENDIF
      ENDIF
!
!     opening the files  *.lsj.lbl
!
      util_csl_file = NAME(1:K-1)//'.lsj'//'.lbl'
      call OPENFL(57,util_csl_file,"formatted","UNKNOWN",ierr)
      IF (IERR .EQ. 1) THEN
         print *, 'Error when opening',util_csl_file
         STOP
      ENDIF
!GG-2017
      IF (UNIQUE .EQ. 1) THEN
         util_csl_file = NAME(1:K-1)//'.uni.lsj'//'.lbl'
         call OPENFL(73,util_csl_file,"formatted","UNKNOWN",ierr)
         IF (IERR .EQ. 1) THEN
            print *, 'Error when opening',util_csl_file
            STOP
         ENDIF
         util_csl_file = NAME(1:K-1)//'.uni.lsj'//'.sum'
         call OPENFL(72,util_csl_file,"formatted","UNKNOWN",ierr)
         IF (IERR .EQ. 1) THEN
            print *, 'Error when opening',util_csl_file
            STOP
         ENDIF
      END IF
!GG-2017 end
      write(57,'(A53)') " Pos   J   Parity      Energy Total      Comp. of ASF"
!
!     Opening the files   *.lsj.j and
!
      IF(ioutj == 1) THEN
         util_csl_file = NAME(1:K-1)//'.lsj'//'.j'
         call OPENFL(58,util_csl_file,"formatted","UNKNOWN",ierr)
         IF (IERR .EQ. 1) THEN
            print *, 'Error when opening',util_csl_file
            STOP
         END IF
         WRITE (58,'(2X,A6,A,F5.1,A,I3,A,I7)' )    &
         NAME(1:K-1),'  Z = ',Z ,'  NEL = ',NELEC,'   NCFG ='! ,asf_set_LS%csf_set_LS%nocsf
      ENDIF
!
!     Opening the files   *.lsj.T and
!
      IF(ioutT == 1) THEN
         util_csl_file = NAME(1:K-1)//'.lsj'//'.T'
!GGR         OPEN(59,FILE=util_csl_file,FORM='unformatted',STATUS='NEW',   &
         OPEN(59,FILE=util_csl_file,FORM='formatted',STATUS='NEW',   &
                                                           IOSTAT=IERR)
         if (ierr /= 0) then
            print *, 'Error when opening ',util_csl_file
            stop
         end if
!GGR         write (59) 'jj2lsj'
         write (59,'(A6)') 'jj2lsj'
!GGR         write (59) NELEC, NCF, NW, NBLOCK
         write (59,'(4I12)') NELEC, NCF, NW, NBLOCK
         DO JB = 1, NBLOCK
             IATJPDUM = TWO_J(JB) +1
!GGR            write (59) JB, NCFINBLK(JB), NEVINBLK(JB), IATJPDUM
            write (59,'(4I12)') JB, NCFINBLK(JB), NEVINBLK(JB), IATJPDUM
         END DO
      ELSE IF(ioutT == 2) THEN
         util_csl_file = NAME(1:K-1)//'.lsj'//'.T'
!GGR         OPEN(59,FILE=util_csl_file,FORM='unformatted',STATUS='OLD',   &
         OPEN(59,FILE=util_csl_file,FORM='formatted',STATUS='OLD',   &
                                                           IOSTAT=IERR)
         IF (IERR /= 0) THEN
            print *, 'Error when opening ',util_csl_file
            stop
         END IF
!GGR         READ (59, IOSTAT=IERR) G92MIX
         READ (59, '(A6)', IOSTAT=IERR) G92MIX
         IF (IERR/=0 .OR. G92MIX/='jj2lsj') THEN
            WRITE (ISTDE, *) 'Not a jj2lsj Transformation File;'
            close(59)
            stop
         ENDIF
!GGR         READ (25) NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
         read (59,'(4I12)') NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
         if(NELECDUM /= NELEC .or. NCFTOTDUM /= NCF .or. NWDUM /= NW   &
                                        .or. NBLOCKDUM /= NBLOCK) then
            print*, NELEC, NCF, NW, NBLOCK
            print*, NELECDUM, NCFTOTDUM, NWDUM, NBLOCKDUM
            print*, "Wrong transformation file *.lsj.T"
            close(59)
            stop
         end if
         DO JB = 1, NBLOCKDUM
            IATJPDUM = TWO_J(JB) +1
!GGR            read (59) JBDUM, NCFINBLKDUM, NEVINBLKDUM, IATJPDUM
            read (59,'(4I12)') JBDUM, NCFINBLKDUM, NEVINBLKDUM, IATJPDUM
            if(JBDUM /= JB .or. NCFINBLKDUM /= NCFINBLK(JB) .or.       &
               NEVINBLKDUM /= NEVINBLK(JB) .or.                        &
                                         IATJPDUM /= TWO_J(JB) +1) then
               print*, JB,NCFINBLK(JB), NEVINBLK(JB), TWO_J(JB) +1
               print*, JBDUM, NCFINBLKDUM, NEVINBLKDUM, IATJPDUM
               print*, "Wrong transformation file *.lsj.T"
               close(59)
               stop
            end if
         END DO
      ENDIF
      END SUBROUTINE inscreen
!
!***********************************************************************
!                                                                      *
      SUBROUTINE inscreenlev(record,levels,number_of_levels,fail)
!                                                                      *
!     Attempts to interprete the serial level numbers which are given  *
!     in record. These level numbers can be given in the format:       *
!                1 3 4  7 - 20  48  69 - 85                            *
!     Any order and 'overlapping' intervals are also supported.        *
!     The procedure returns with fail = .true. if the level numbers    *
!     cannot be interpreted properly (fail = .false. otherwise).       *
!                                                                      *
!     The level numbers are returned in the vector                     *
!                                          levels(1:number_of_levels). *
!     The procedure assumes that this vector has a sufficient          *
!     dimension to store all level numbers.                            *
!                                                                      *
!     Calls: idigit.                                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE idigit_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*), intent(in) :: record
      logical, intent(out)         :: fail
      integer, intent(out)         :: number_of_levels
      integer, dimension(1:Vectors_number), intent(inout) :: levels
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!GG      logical, dimension(200)   :: low_to
      logical, dimension(2000)   :: low_to
      character(len=500)        :: string
      integer                   :: a, i, lower, n
      integer, dimension(2000)   :: low
!GG      integer, dimension(200)   :: low
      integer, dimension(Vectors_number) :: run
!-----------------------------------------------
      fail = .true.;  levels(1:Vectors_number) = 0;  number_of_levels = 0
      run(1:Vectors_number) = 0;   string = adjustl(record)
      n = 0;   lower = 0;   low_to(:) = .false.
    1 string = adjustl(string)
      if (string(1:1) == " ") then
         if (n == 0  .or.  low_to(n)) then
            return
         else
            goto 10
         end if
      else if (string(1:1) == "-") then
         if (n == 0)  return
         low_to(n)   = .true.;    string(1:1) = " "
         goto 1
      else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or.  &
               string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
               string(1:1) == "4"   .or.   string(1:1) == "5" .or.  &
               string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
               string(1:1) == "8"   .or.   string(1:1) == "9") then
         a = idigit(string(1:1));    lower = 10*lower + a
         if (string(2:2) == " "   .or.   string(2:2) == "-") then
            n = n + 1;    low(n) = lower;    lower = 0
         end if
         string(1:1) = " "
         goto 1
      end  if
!
!     Determine no_eigenpairs and max_eigenpair
   10 levels(:) = 0
      do  i = 1,n
         if (low_to(i)) then
            if (low(i) <= low(i+1)) then;   run(low(i):low(i+1)) = 1
            else;                           run(low(i+1):low(i)) = 1
            end if
            if (low_to(i+1)) then;   return
            else;                    cycle
            end if
         else
            run(low(i)) = 1
         end if
      end do
      number_of_levels = 0
      do  i = 1,Vectors_number
         if (run(i) == 1) then
            number_of_levels         = number_of_levels + 1
            levels(number_of_levels) = i
         end if
      end do
      fail = .false.
      END SUBROUTINE inscreenlev
!
!***********************************************************************
!                                                                      *
      SUBROUTINE jj2lsj
!                                                                      *
!     Controls the transformation of atomic states from a jj-          *
!     to a LS-coupled CSF basis.                                       *
!                                                                      *
!     Calls: asf2LS, convrt_double, dallocASFLS, inscreen, ispar,      *
!            itjpo, packlsCSF, prCSFLSall, prCSFjj, prCSFLS, prCSFall, *
!            setLS.                                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!     Modified by G. Gaigalas,                              May 2021   *
!     Modified by G. Gaigalas                                   2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE EIGV_C,          ONLY: EAV, EVAL, EVEC
      USE ORB_C,           ONLY: NCF
      USE PRNT_C,          ONLY: IVEC, NVEC
      USE IOUNIT_C,        ONLY: ISTDE
      USE CONS_C,          ONLY: ZERO, ONE
      USE def_C,           ONLY: Z, NELEC
      USE blk_C,           ONLY: NEVINBLK, NCFINBLK, NBLOCK, TWO_J
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE ispar_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!GG NIST
      integer :: i, j, jj, ii, string_l,IBLKNUM,ioutC,ioutj,UNIQUE,ioutT
      integer :: level, nocsf_min, lev, string_length
      integer :: nocsf_max, sum_nocsf_min, Before_J
!GG NIST
      integer :: LOC, NCFMIN, NCFMAX, NCF_LS_jj_MAX
      real(DOUBLE) :: THRESH, wa, wb,Ssms,g_j,g_JLS, sumthrsh
      character(len=4)   :: string_CNUM
!GG      character(len=64)  :: string_CSF_ONE
      character(len=164)  :: string_CSF_ONE
      integer,  dimension(Blocks_number)      :: number_of_levels
      integer, dimension(Blocks_number,Vectors_number) :: levels
      integer, dimension(1:Vectors_number)     :: max_comp
!GG      integer, dimension(1:Vectors_number)     :: leading_LS
      integer, dimension(:), pointer     :: leading_LS
!GG      integer, dimension(1:100)       :: iw
      integer, dimension(:), pointer       :: iw
!GG      real(DOUBLE), dimension(1:100)  :: weights, weights2
      real(DOUBLE), dimension(:), pointer  :: weights, weights2
      integer, dimension(:), pointer  :: ithresh
!GG      character(LEN=64), dimension(1:Vectors_number) :: string_CSF
      character(LEN=164), dimension(1:Vectors_number) :: string_CSF
!-----------------------------------------------
      Ssms = ZERO;   g_j = ZERO;   g_JLS = ZERO;    Before_J = 0
      call inscreen(THRESH,levels,number_of_levels,ioutC,ioutj,UNIQUE, &
                                                                 ioutT)
      allocate(ithresh(NCF))
      do  IBLKNUM = 1, NBLOCK
         if(IBLKNUM == 1) THEN
            NCFMIN = 1
            NCFMAX = NCFINBLK(IBLKNUM)
         else
            NCFMIN = NCFMAX + 1
            NCFMAX = NCFMIN + NCFINBLK(IBLKNUM) - 1
         end if
         if(number_of_levels(IBLKNUM) == 0) GO TO 1
         ithresh = 0
         do  i = NCFMIN, NCFMAX
            sumthrsh = ZERO
            do  lev = 1, number_of_levels(IBLKNUM)
                level = levels(IBLKNUM,lev)
                LOC = (level-1)*NCF
                sumthrsh = sumthrsh + dabs(EVEC(i+LOC))
            end do
            if(dabs(sumthrsh) >= dabs(EPSNEW)) ithresh(i) = 1
         end do
         call setLS(ithresh,NCFMIN,NCFMAX)
!
!     output to *.lsj.c
         if(ioutC == 1) call prCSFLSall (56,asf_set_LS%csf_set_LS,Before_J)
         allocate(asf_set_LS%asf(1:NVEC))
         do i = 1, NVEC
            allocate(asf_set_LS%asf(i)%eigenvector(1:asf_set_LS%csf_set_LS%nocsf))
         end do
         do lev = 1, number_of_levels(IBLKNUM)
            level = levels(IBLKNUM,lev)
            asf_set_LS%asf(level)%level_No = level
            asf_set_LS%asf(level)%energy = EAV+EVAL(level)
            asf_set_LS%asf(level)%eigenvector = ZERO
         end do
         asf_set_LS%noasf = number_of_levels(IBLKNUM)
         nocsf_max = asf_set_LS%csf_set_LS%nocsf
!GG NIST
         NCF_LS_jj_MAX = MAX(1000,NCFMAX-NCFMIN+1,asf_set_LS%csf_set_LS%nocsf+1)
         allocate (weights(NCF_LS_jj_MAX))
         allocate (weights2(NCF_LS_jj_MAX))
         allocate (iw(NCF_LS_jj_MAX))
         allocate (leading_LS(NCF_LS_jj_MAX))
         do  lev = 1, number_of_levels(IBLKNUM)
            level = levels(IBLKNUM,lev)
            weights = ZERO;    iw = 0;    wb = ZERO
            do  i = NCFMIN, NCFMAX
               if(ithresh(i) == 1) then
                  LOC = (level-1)*NCF
                  wa = EVEC(i+LOC) * EVEC(i+LOC)
                  wb = wb + wa
                  do  j = 1,999
                     if (wa > weights(j)) then
                        weights(j+1:1000) = weights(j:999)
                        weights(j)       = wa
                        iw(j+1:1000)      = iw(j:999)
                        iw(j)            = i - NCFMIN + 1
                        exit
                     end if
                  end do
               end if
            end do
            if(lev == 1) then
               print *, " "
               WRITE(*,'(A,A)') " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .",&
               "  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
               WRITE(*, '(A,2X,I4,16X,A,I3)')                            &
               " Under investigation is the block:",IBLKNUM,             &
               " The number of eigenvectors:", number_of_levels(IBLKNUM)
               WRITE(*,'(A,I10,10X,A,I10)')                              &
               " The number of CSF (in jj-coupling):",NCFINBLK(IBLKNUM), &
               " The number of CSF (in LS-coupling):",asf_set_LS%csf_set_LS%nocsf
            else
               print *, " "
               WRITE(*,'(A)') " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
               WRITE(*,'(A)') " The new level is under investigation."
            end if
!
!     perform the transformation
!
            if(lev == 1) call asf2ls                                    &
            (iw(1),ithresh,number_of_levels(IBLKNUM),IBLKNUM,levels,   &
             NCFMIN,NCFMAX,ioutT)
!
!     output to the screen jj- coupling
            print *, "Weights of major contributors to ASF in jj-coupling:"
            print *, " "
            print *, " Level  J Parity      CSF contributions"
            print *, " "
            if (level > 1   .and.   dabs(wb) > 1.0001) then
               print *, "level, wb = ",level,wb
               stop "JJ2LSJ(): program stop A."
            end if
            nocsf_min = 5
            do j = 1,5
               if(abs(weights(j)) < 0.00001) then
                  nocsf_min = j - 1
                  exit
               end if
            end do
            IF(ISPAR(iw(1)) == -1) THEN
              asf_set_LS%asf(level)%parity = "-"
            ELSE IF(ISPAR(iw(1)) == 1) THEN
              asf_set_LS%asf(level)%parity = "+"
            ELSE
              STOP "JJ2LSJ: program stop D."
            END IF
            asf_set_LS%asf(level)%totalJ = ITJPO(iw(1)+NCFMIN-1) - 1
            call convrt_double(asf_set_LS%asf(level)%totalJ,string_CNUM,string_l)
            print 16, IVEC(level),string_CNUM(1:string_l),               &
              asf_set_LS%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
            print*, "              Total sum over  weight (in jj) is:",wb
            print *, " "
            print *, "Definition of leading CSF:"
            print *, " "
!CGG            call prCSFjj(-1,iw(1))
            call prCSFjj(-1,iw(1)+NCFMIN-1,iw(1))
!
!     output to the screen LS- coupling
            print *, " "
            print *, " "
            print *, "Weights of major contributors to ASF in LS-coupling:"
            print *, " "
            print *, " Level  J Parity      CSF contributions"
            print *, " "
            sum_nocsf_min = 0
            weights = ZERO;  weights2 = ZERO;  iw = 0;  wb = ZERO
            do  i = 1,asf_set_LS%csf_set_LS%nocsf
               wa = asf_set_LS%asf(level)%eigenvector(i)
               wb = wb + wa*wa
!GG NIST
               if(i == 1) then
                  weights(1)  = wa
                  weights2(1) = wa*wa
                  iw(1)       = 1
               else
                  do  j = 1,i
                     if(j == i) then
                       weights(i)  = wa
                       weights2(i) = wa*wa
                       iw(i)       = i
                     else
                       if (wa*wa > weights2(j)) then
                         do ii = i,j,-1
                           weights2(ii+1) = weights2(ii)
                           weights(ii+1)  = weights(ii)
                           iw(ii+1)       = iw(ii)
                         end do
                         weights(j)  = wa
                         weights2(j) = wa*wa
                         iw(j)       = i
                         exit
                       end if
                     end if
                  end do
               end if
            end do
            if(IMINCOMPOFF == 1) THEN
               if (level > 1 .and. (wb*100.+dabs(MINCOMP)) <= 100.) THEN
                  WRITE (*,"(A,F5.0,A)")   &
                    "Attention!!!  The program not reach the accuracy",MINCOMP," %"
                  print *, "level, sum of weights in LSJ = ",level,wb
                  print *,""
               end if
            else
               if (level > 1 .and. abs(wb) > 1.001) then
                  print*,"Attention!!!  The sum of weights is bigger then 1.0001"
                  print *, "level, sum of weights in LSJ = ",level,wb
               end if
               if (level > 1 .and.  abs(wb) < .998) then
                  print*,"Attention!!!  The sum of weights is little then .998"
                  print *, "level, sum of weights in LSJ = ",level,wb
               end if
            end if
            nocsf_min = nocsf_max;   jj = 0
            do j = 1,nocsf_max
               if(dabs(weights2(j)) < dabs(THRESH)) then
                  nocsf_min = j - 1
                  exit
               end if
               jj = jj + 1
               leading_LS(jj) = iw(j)
            end do
            max_comp(level) = IW(1)
            sum_nocsf_min = sum_nocsf_min + nocsf_min
            if (nocsf_min <= 0) then
               stop "JJ2LSJ(): program stop C."
            end if
            call convrt_double(asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,string_CNUM,string_l)
            print 16, asf_set_LS%asf(level)%level_No,string_CNUM(1:string_l), &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,                   &
                   (weights2(j),iw(j),J=1,5)
            print*, "              Total sum over  weight (in LSJ) is:",wb
            print *, " "
            print *, "Definition of leading CSF:"
            print *, " "
            if(sum_nocsf_min > 5) sum_nocsf_min = 5
            do  i = 1,NCF
               do  j = 1, sum_nocsf_min
                  if (i == iw(j)) then
                     call prCSFLS (-1,asf_set_LS%csf_set_LS,leading_LS(j))
                     exit
                  end if
               end do
            end do
!
!     output to *.lsj.lbl
            IF(Before_J == 0) THEN
               Before_J = 1
            ELSE IF(lev == 1 .and. Before_J == 1) THEN
               write(57,'(A1)') ' '
            END IF
            DO j = 1,nocsf_min
               CALL packlsCSF(asf_set_LS%csf_set_LS,iw(j),string_CSF_ONE)
               IF(J == 1) THEN
                  write(57,'(1X,I2,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)')     &
                  asf_set_LS%asf(level)%level_No,string_CNUM(1:string_l),&
                  asf_set_LS%csf_set_LS%csf(iw(1))%parity,               &
                  asf_set_LS%asf(level)%energy,wb*100,"%"
                  string_CSF(level) = string_CSF_ONE
               END IF
               string_length = Len_Trim(string_CSF_ONE)
               write(57,'(7X,F12.8,3X,F11.8,3X,A)') weights(j),weights2(j),string_CSF_ONE(1:string_length)
            END DO
!     output  *.lsj.j and
            IF(ioutj == 1) THEN
               IF(lev == 1) THEN
                  WRITE (58, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',       &
                  asf_set_LS%asf(level)%totalJ,'NUMBER =',number_of_levels(IBLKNUM)
               END IF
               string_length = Len_Trim(string_CSF(level))
               write(58,'(3(A8,F15.10))') 'Ssms=', Ssms, 'g_J=',g_J, 'g_JLS=', g_JLS
               write(58,'(I6,F16.9,2X,A)')                                &
                 max_comp(level),asf_set_LS%asf(level)%energy,string_CSF(level)(1:string_length)
               write(58,'(7F11.8)')                                       &
                 (asf_set_LS%asf(level)%eigenvector(i), i=1,asf_set_LS%csf_set_LS%nocsf)
            END IF
         END DO
         call dallocASFLS(asf_set_LS)
         deallocate(weights)
         deallocate(weights2)
         deallocate(iw)
         deallocate(leading_LS)
    1    CONTINUE
      END DO
      deallocate(ithresh)
      close(56)
!GG-2017
      IF (UNIQUE .EQ. 1) THEN
         CALL uniquelsj
         close(72)
         close(73)
      end if
!GG-2017
      close(57)
      IF(ioutj == 1) write(58,'(A)') '**'
      close(58)
      deallocate(Jcoup)
   16 format(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5,&
             " of",i7,3x,f8.5," of",i7/16X))
      return
      END SUBROUTINE jj2lsj
!
!***********************************************************************
!                                                                      *
      SUBROUTINE packlsCSF(csf_set_LS,csf_number,string_CSF)
!                                                                      *
!     Encording all CSF in LS - coupling.                              *
!                                                                      *
!     Rules for encoding                                               *
!        1. All blanks deleted                                         *
!        2. If Qi=1, omit Qi                                           *
!        3. If Qi=1 or Qi>=4l+1, omit ALFAi                            *
!        4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.     *
!                                                                      *
!     Calls: getchLS, packls.                                          *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      type(csf_basis_LS), intent(in) :: csf_set_LS
      integer, intent(in)            :: csf_number
!GG      CHARACTER(LEN=64), INTENT(OUT) :: string_CSF
      CHARACTER(LEN=164), INTENT(OUT) :: string_CSF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: I, j, counter
      integer, dimension(:), pointer          :: occupation
!GG      integer, dimension(8)                   :: Q
      integer, dimension(18)                   :: Q
      character(len=4)                        :: LS, XLS
      character(len=4), dimension(:),pointer  :: string_LS, string_XLS
!GG      character(len=3), dimension(15)         :: COUPLE
      character(len=3), dimension(35)         :: COUPLE
!GG      character(len=3), dimension(8)          :: ELC
      character(len=3), dimension(18)          :: ELC
      CHARACTER(LEN=2)                        :: String
!GG      CHARACTER(LEN=64)                       :: string_CSF_PACK
      CHARACTER(LEN=164)                       :: string_CSF_PACK
      CHARACTER(LEN=3), DIMENSION(:), POINTER :: String_NL
      CHARACTER(LEN=1), DIMENSION(0:20)       :: L1
      DATA L1 /'s','p','d','f','g','h','i','k','l','m','n', &
               'o','q','r','t','u','v','w','x','y','z'/
!-----------------------------------------------
      counter = 0
      allocate(occupation(csf_set_LS%nwshells))
      allocate(string_LS(csf_set_LS%nwshells))
      allocate(string_XLS(csf_set_LS%nwshells))
      do j = csf_set_LS%nwcore+1, csf_set_LS%nwshells
         if (csf_set_LS%csf(csf_number)%occupation(j) > 0) then
            counter = counter +1;        occupation(counter) = j
            call getchLS(csf_number,j,LS,XLS)
            string_LS(counter)  = LS;    string_XLS(counter) = XLS
         end if
      end do
!
      allocate(String_NL(counter))
      DO I = 1,counter
         write(String,'(I2)')csf_set_LS%shell(occupation(I))%n
         String_NL(I)(1:2) = String(1:2)
         String_NL(I)(3:3) = L1(csf_set_LS%shell(occupation(I))%l)
      END DO
!GG      IF(counter > 8) THEN
      IF(counter > 18) THEN
            print*, counter
            STOP "prCSFLS: A"
      END IF
      DO I = 1,counter
         COUPLE(I)(1:3) = string_LS(I)(1:3)
         IF(I > 1) THEN
            COUPLE(I+counter-1)(1:3) = string_XLS(I)(1:3)
         END IF
         Q(I) = csf_set_LS%csf(csf_number)%occupation(occupation(I))
         ELC(I)(1:3) = String_NL(I)(1:3)
      END DO
      CALL PACKLS(counter,ELC,Q,COUPLE,string_CSF_pack)
      string_CSF = string_CSF_pack
      deallocate(String_NL)
      deallocate(occupation)
      deallocate(string_LS)
      deallocate(string_XLS)
      END SUBROUTINE packlsCSF
!
!***********************************************************************
!                                                                      *
      SUBROUTINE prCSFjj(stream,csf_number,csf_number_local)
!                                                                      *
!     Print all information about the single CSF scheme csf_set in     *
!     a nead format on stream.                                         *
!                                                                      *
!     Calls: convrt_double, itjpo, iq, jqs.                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE ORB_C,           ONLY: NW, NP, NAK
      USE M_C,             ONLY: NCORE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE IQ_I
      USE JQS_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: stream, csf_number, csf_number_local
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: I, j, counter, string_lenth
      integer, dimension(:), pointer          :: occupation
      character(len=4), dimension(:), pointer :: string_J, string_X
      character(len=4)                        :: string_CNUM
      CHARACTER(LEN=2)                        :: String
      CHARACTER(LEN=4), DIMENSION(:), POINTER ::  String_NL
      CHARACTER(LEN=1), DIMENSION(0:20)       :: L1
      DATA L1 /'s','p','d','f','g','h','i','k','l','m','n', &
               'o','q','r','t','u','v','w','x','y','z'/
!-----------------------------------------------
      counter = 0
      allocate(occupation(NW))
      allocate(string_J(NW))
      allocate(string_X(NW))
      do j = NCORE+1, NW
         if (IQ(j,csf_number) > 0) then
            counter = counter +1
            occupation(counter) = j
            if (JQS(3,j,csf_number)-1 == 0) then
               string_J(counter) = "    "
            else
               call convrt_double(JQS(3,j,csf_number)-1,string_CNUM,string_lenth)
               string_J(counter) = string_CNUM(1:string_lenth)
            end if
            if (Jcoup(j,csf_number) == 0) then
               string_X(counter) = "     "
            else
               call convrt_double(Jcoup(j,csf_number),string_CNUM,string_lenth)
               string_X(counter) =  string_CNUM(1:string_lenth)
            end if
         end if
      end do
!
      allocate(String_NL(counter))
      DO I = 1,counter
         write(String,'(I2)')NP(occupation(I))
         String_NL(I)(1:2) = String(1:2)
         J = ((IABS(NAK(occupation(I)))*2)-1+    &
             NAK(occupation(I))/IABS(NAK(occupation(I))))/2
         String_NL(I)(3:3) = L1(J)
         IF(NAK(occupation(I)) <= 0) THEN
         String_NL(I)(4:4) = " "
         ELSE
         String_NL(I)(4:4) = "-"
         END IF
      END DO
!
      if (stream == -1) then
         write(*,1) csf_number_local,(String_NL(j),                   &
            IQ(occupation(j),csf_number),j=1,counter)
         write(*,2)(string_J(j),j=1,counter)
         call convrt_double(ITJPO(csf_number)-1,string_CNUM,string_lenth)
         write(*,3)(string_X(j),j=2,counter-1),string_CNUM(1:string_lenth)
      end if
      deallocate(occupation)
      deallocate(string_J)
      deallocate(string_X)
    1 format(4x,i9,')',100(a4,'(',i2,')',2x))
    2 format(18x,a4,100(6X,a4))
    3 format(31x,a4,100(6X,a5))
      END SUBROUTINE prCSFjj
!
!***********************************************************************
!                                                                      *
      SUBROUTINE prCSFLS(stream,csf_set_LS,csf_number)
!                                                                      *
!     Print all information about the single CSF scheme csf_set in     *
!     a nead format on stream.                                         *
!                                                                      *
!     Calls: convrt_double, getchLS.                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)            :: stream, csf_number
      type(csf_basis_LS), intent(in) :: csf_set_LS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: I, j, counter, string_lenth
      integer, dimension(:), pointer          :: occupation
      character(len=4)                        :: LS, XLS
      character(len=4), dimension(:),pointer  :: string_LS, string_XLS
      character(len=4)                        :: string_CNUM
      CHARACTER(LEN=2)                        :: String
      CHARACTER(LEN=3), DIMENSION(:), POINTER :: String_NL
      CHARACTER(LEN=1), DIMENSION(0:20)       :: L1
      DATA L1 /'s','p','d','f','g','h','i','k','l','m','n', &
               'o','q','r','t','u','v','w','x','y','z'/
!-----------------------------------------------
      counter = 0
      allocate(occupation(csf_set_LS%nwshells))
      allocate(string_LS(csf_set_LS%nwshells))
      allocate(string_XLS(csf_set_LS%nwshells))
      do j = csf_set_LS%nwcore+1, csf_set_LS%nwshells
         if (csf_set_LS%csf(csf_number)%occupation(j) > 0) then
            counter = counter +1;        occupation(counter) = j
            call getchLS(csf_number,j,LS,XLS)
            string_LS(counter)  = LS;    string_XLS(counter) = XLS
         end if
      end do
!
      allocate(String_NL(counter))
      DO I = 1,counter
         write(String,'(I2)')csf_set_LS%shell(occupation(I))%n
         String_NL(I)(1:2) = String(1:2)
         String_NL(I)(3:3) = L1(csf_set_LS%shell(occupation(I))%l)
      END DO
!
      if (stream == -1) then
         write(*,3)csf_number,(String_NL(j),                                 &
            csf_set_LS%csf(csf_number)%occupation(occupation(j)),j=1,counter)
         call convrt_double(csf_set_LS%csf(csf_number)%totalJ,string_CNUM,string_lenth)
         write(*,4)(string_LS(j),j=1,counter),(string_XLS(j),j=2,counter),&
              string_CNUM(1:string_lenth)
      else
         write(stream,1)(String_NL(j),                         &
            csf_set_LS%csf(csf_number)%occupation(occupation(j)),j=1,counter)
         call convrt_double(csf_set_LS%csf(csf_number)%totalJ,string_CNUM,string_lenth)
         write(stream,2)(string_LS(j),j=1,counter),(string_XLS(j),j=2,counter)
      end if
      deallocate(String_NL)
      deallocate(occupation)
      deallocate(string_LS)
      deallocate(string_XLS)
    1 format(8(1X,A3,'(',i2,')'))
    2 format(1X,15(A4))
    3 format(i10,')',4x,100(a3,'(',i2,')',2x))
    4 format(19x,a4,100(5X,a4))
      END SUBROUTINE prCSFLS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE prCSFLSall(stream,csf_set_LS,Before_J)
!                                                                      *
!     Print all information about the CSF scheme csf_set in a nead     *
!     format on stream.                                                *
!                                                                      *
!     Calls: prCSFLS.                                                  *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)            :: stream
      type(csf_basis_LS), intent(in) :: csf_set_LS
      integer, intent(in)            :: Before_J
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                  :: I, J
      CHARACTER (LEN=2)                        :: String
      CHARACTER (LEN=3), DIMENSION(:), POINTER :: String_NL
      CHARACTER (LEN=1), DIMENSION(0:20)       :: L1
      DATA L1 /'s','p','d','f','g','h','i','k','l','m','n', &
               'o','q','r','t','u','v','w','x','y','z'/
!-----------------------------------------------
      if(Before_J == 0) then
         write(stream,*) " "
         allocate(String_NL(csf_set_LS%nwcore))
         do I = 1,csf_set_LS%nwcore
            write(String,'(I2)')csf_set_LS%shell(I)%n
            String_NL(I)(1:2) = String(1:2)
            String_NL(I)(3:3) = L1(csf_set_LS%shell(I)%l)
         end do
         write(stream,1)(String_NL(J),J=1,csf_set_LS%nwcore)
         deallocate(String_NL)
      end if
!
!     output to *.lsj.c
      do i = 1,csf_set_LS%nocsf
         call prCSFLS(stream,csf_set_LS,i)
      end do
      write(stream,*) "*"
    1 format(18(1X,A3))
      END SUBROUTINE prCSFLSall
!
!***********************************************************************
!                                                                      *
      SUBROUTINE setLS(ithresh,NCFMIN,NCFMAX)
!                                                                      *
!     This subroutine fills up the variable asf_set_LS%csf_set_LS      *
!     with data generated using the one from asf_set%csf_set.          *
!                                                                      *
!     This subroutine contains the following internal subroutines:     *
!          * subroutine            setLS_action                        *
!          * subroutine            setLS_add_quantum_numbers           *
!          * recursive subroutine  setLS_job_count                     *
!          * function              setLS_equivalent_csfs               *
!                                                                      *
!     Calls: itjpo, iq, setLS_equivalent_csfs, setLS_job_count.        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE DEF_C,           ONLY: NELEC
      USE ORB_C,           ONLY: NW, NP, NAK, NCF
      USE m_C,             ONLY: NCORE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE IQ_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, dimension(:), intent(in)  :: ithresh
      integer, intent(in)  :: NCFMIN, NCFMAX
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(nl), dimension(:), pointer :: shell_temp
      integer, dimension(:), pointer  :: nlLSval, nljjval
      type(lsj_list)                  :: nonequiv_csfs_jj
!
      integer :: isubc, isubc2, icsf_jj, icsf_jj2, icsf_jj_real, icsf_LS
      integer :: N, action_type
      logical :: new_one, found_parent_minus, found_parent_plus
!
      integer, dimension(:), pointer :: all_occupation
      integer, dimension(:), pointer :: Li, Si, L_i, S_i, w, Q
      integer                        :: J
!-----------------------------------------------
      asf_set_LS%csf_set_LS%number_of_electrons = NELEC
!
!     1. define  nl, parent
      allocate(shell_temp(NW));  allocate(nlLSval(NW));  allocate(nljjval(NW))
      nljjval = 0;  nlLSval = 0
!
      asf_set_LS%csf_set_LS%nwshells = 0
      do isubc = 1, NW, 1
         nljjval(isubc) = NP(isubc)*100 +                              &
                  ((IABS(NAK(isubc))*2)-1+NAK(isubc)/IABS(NAK(isubc)))/2
         new_one = .true.
         do isubc2 = 1, isubc-1, 1
            if(nljjval(isubc) == nljjval(isubc2)) new_one = .false.
         end do
         if(new_one) then
            asf_set_LS%csf_set_LS%nwshells = asf_set_LS%csf_set_LS%nwshells + 1
            shell_temp(asf_set_LS%csf_set_LS%nwshells)%n = NP(isubc)
            shell_temp(asf_set_LS%csf_set_LS%nwshells)%l=                &
            ((IABS(NAK(isubc))*2)-1+NAK(isubc)/IABS(NAK(isubc)))/2
            nlLSval(asf_set_LS%csf_set_LS%nwshells) = nljjval(isubc)
         end if
      end do
!
      allocate(asf_set_LS%csf_set_LS%shell(asf_set_LS%csf_set_LS%nwshells))
      allocate(asf_set_LS%csf_set_LS%parent(asf_set_LS%csf_set_LS%nwshells))
!
      do isubc=1, asf_set_LS%csf_set_LS%nwshells, 1
         asf_set_LS%csf_set_LS%shell(isubc) = shell_temp(isubc)
      end do
!
      deallocate(shell_temp)
!
!     begin find parent
      do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
         found_parent_minus = .false.;  found_parent_plus = .false.
         do isubc2 = 1, NW, 1
            if(nlLSval(isubc) == nljjval(isubc2)) then
               if(NAK(isubc2) > 0) then
                  found_parent_minus = .true.
                  asf_set_LS%csf_set_LS%parent(isubc)%parent_minus = isubc2
               else
                  found_parent_plus = .true.
                  asf_set_LS%csf_set_LS%parent(isubc)%parent_plus  = isubc2
               end if
            end if
         end do
         if(.not.found_parent_plus)  &
            asf_set_LS%csf_set_LS%parent(isubc)%parent_plus = 0
         if(.not.found_parent_minus)  &
            asf_set_LS%csf_set_LS%parent(isubc)%parent_minus = 0
      end do
!
!     end find parent  --------------------
!
!     2. define the number of "core" shells
!        (the LS shell is supposed to be "core" if:
!        1. l=0 and corresponding jj subshell is "core"
!        2. l<>0 l+ and l - "core" subshells            )
      asf_set_LS%csf_set_LS%nwcore=0
      do isubc=1, asf_set_LS%csf_set_LS%nwshells, 1
         if(asf_set_LS%csf_set_LS%parent(isubc)%parent_minus .le. NCORE &
            .and.                                                       &
            asf_set_LS%csf_set_LS%parent(isubc)%parent_plus .le. NCORE) &
            asf_set_LS%csf_set_LS%nwcore = asf_set_LS%csf_set_LS%nwcore+1
      end do
!
!     3. form the list of "nonequivalent" csfs_jj
!        (i.e. csfs_jj different in J,parity, or
!         some l's occupation numbers Ni = N_(i+) + N_(i-))
      allocate(nonequiv_csfs_jj%items(NCF))
      nonequiv_csfs_jj%list_size = 0
      do icsf_jj = NCFMIN, NCFMAX
         new_one = .true.
         if(ithresh(icsf_jj) == 1) then
            do icsf_jj2 = NCFMIN , icsf_jj - 1
               if(ithresh(icsf_jj2) == 1) then
                  if(setLS_equivalent_csfs(icsf_jj,icsf_jj2)) then
                     new_one = .false.
                     exit
                  end if
               end if
            end do
            if(new_one) then
               nonequiv_csfs_jj%list_size = nonequiv_csfs_jj%list_size + 1
               nonequiv_csfs_jj%items(nonequiv_csfs_jj%list_size) = icsf_jj
            end if
         end if
      end do
!
!     4. for each nonequivalent csf_jj find all the csfs_LS
!
!        To avoid the dependency on the number of subshells
!         the recursive subroutine is used
      allocate(Li(asf_set_LS%csf_set_LS%nwshells))
      allocate(L_i(asf_set_LS%csf_set_LS%nwshells))
      allocate(Si(asf_set_LS%csf_set_LS%nwshells))
      allocate(S_i(asf_set_LS%csf_set_LS%nwshells))
      allocate(Q(asf_set_LS%csf_set_LS%nwshells))
      allocate(w(asf_set_LS%csf_set_LS%nwshells))
      allocate(all_occupation(asf_set_LS%csf_set_LS%nwshells))
!
!     4.1 - find the number of csfs_LS
      asf_set_LS%csf_set_LS%nocsf = 0
!
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
!
!     set the variable for conviniency ...
         icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         J = ITJPO(icsf_jj_real)-1
!
!     define the occupation numbers
!
         do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
            all_occupation(isubc) = 0
            do isubc2 = 1, NW, 1
               if(nlLSval(isubc) == nljjval(isubc2))       &
               all_occupation(isubc)=all_occupation(isubc)+IQ(isubc2,icsf_jj_real)
            end do  !isubc2
!
         end do  !isubc
         N = 0;    action_type = 1
         call setLS_job_count(1, N)
         asf_set_LS%csf_set_LS%nocsf=asf_set_LS%csf_set_LS%nocsf + N
      end do
!
!     4.2 - fill up the csf_LS arrays with the corresponding quantum numbers
!     4.2.1  allocate the arrays
      allocate(asf_set_LS%csf_set_LS%csf(asf_set_LS%csf_set_LS%nocsf))
      do N = 1, asf_set_LS%csf_set_LS%nocsf, 1
         allocate(asf_set_LS%csf_set_LS%csf(N)%occupation(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%seniority(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%w(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellL(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellS(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellLX(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellSX(asf_set_LS%csf_set_LS%nwshells))
      end do
!
!     4.2.2 - fill them with quantum numbers
      icsf_LS = 0
!
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
!
!     set the variable for conveniency ...
         icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         J = ITJPO(icsf_jj_real)-1
         do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
            all_occupation(isubc)=0
            do isubc2 = 1, NW, 1
               if(nlLSval(isubc) == nljjval(isubc2))            &
               all_occupation(isubc)=all_occupation(isubc)+IQ(isubc2,icsf_jj_real)
            end do  !isubc2
         end do  !isubc
         action_type = 2
         call setLS_job_count(1, N)
      end do
!
      deallocate(nonequiv_csfs_jj%items)
      deallocate(nlLSval);  deallocate(nljjval)
!
      deallocate(Q);   deallocate(w);    deallocate(all_occupation)
      deallocate(Li);  deallocate(L_i);  deallocate(Si);
      deallocate(S_i)
!
      CONTAINS
!
!***********************************************************************
!                                                                      *
      SUBROUTINE setLS_action(tip,irez)
!                                                                      *
!     The subroutine defines the "action" of subroutine                *
!     setLS_job_count: whether it counts the number                    *
!     of csfs_LS (asf_set_LS%csf_set_LS%novcsf) or fills the           *
!     arrays of wave functions in LS coupling with                     *
!     asf_set_LS%csf_set_LS%csf(...) with the corresponding            *
!     quantum nubers.                                                  *
!                                                                      *
!     Calls: setLS_add_quantum_numbers.                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: tip, irez
!-----------------------------------------------
      if (tip .eq. 1) then
         irez = irez + 1
      else
         call setLS_add_quantum_numbers()
      end if
      END SUBROUTINE setLS_action
!
!***********************************************************************
!                                                                      *
      SUBROUTINE setLS_add_quantum_numbers()
!                                                                      *
!     The subroutine adds quantum numbers stored in temprorary arrays  *
!     Li, Si, L_i, S_i, w, Q to the corresponding arrays               *
!     of asf_set_LS%csf_set_LS%csf() (i.e. to the arrays of            *
!     corresponding quantum numbers of the wave function in LS         *
!     coupling).                                                       *
!                                                                      *
!     Calls: gettermLS, ispar setLS_action, setLS_job_count.           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ispar_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isubcx
!-----------------------------------------------
      icsf_LS = icsf_LS + 1
!
      if(icsf_LS  .gt.  asf_set_LS%csf_set_LS%nocsf) then
         stop 'setLS_add_quantum_numbers(): program stop A.'
      end if
!
      asf_set_LS%csf_set_LS%csf(icsf_LS)%totalJ = J
!
      IF (ISPAR(icsf_jj_real) == 1) THEN
         asf_set_LS%csf_set_LS%csf(icsf_LS)%parity = "+"
      ELSE IF (ISPAR(icsf_jj_real) == -1) THEN
         asf_set_LS%csf_set_LS%csf(icsf_LS)%parity = "-"
      END IF
!
      do isubcx = 1, asf_set_LS%csf_set_LS%nwshells, 1
         asf_set_LS%csf_set_LS%csf(icsf_LS)%occupation(isubcx) =        &
                                             all_occupation(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%shellL(isubcx)  = Li(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%shellS(isubcx)  = Si(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%shellLX(isubcx) = L_i(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%shellSX(isubcx) = S_i(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%w(isubcx)       = w(isubcx)
         asf_set_LS%csf_set_LS%csf(icsf_LS)%seniority(isubcx) =         &
              2*asf_set_LS%csf_set_LS%shell(isubcx)%l + 1 - Q(isubcx)
      end do
      END SUBROUTINE setLS_add_quantum_numbers
!
!***********************************************************************
!                                                                      *
      RECURSIVE SUBROUTINE setLS_job_count(isubc, rez)
!                                                                      *
!     Recursive subroutine for the calculation of the                  *
!     number of csfs_LS and corresponding quantum numbers.             *
!                                                                      *
!     Calls: setLS_action, setLS_job_count, gettermLS, ittk.           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE jj2lsj_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: isubc, rez
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(subshell_term_LS), dimension(120) :: LS_terms
      integer :: iterm, nr_terms, suma, delta_J
!GG 2015.12.02 NIST
!GG      integer :: number
      integer :: number, numberGG
      integer :: tempLmax,tempLmin,tempSmax,tempSmin,tempS,tempL
      integer :: l_shell, N
!-----------------------------------------------
      if(isubc.gt.(asf_set_LS%csf_set_LS%nwshells) .or. isubc.lt.1) then
         print *, 'isubc = ', isubc
         stop "setLS_job_count(): program stop A."
      end if
!
      if(isubc.le.asf_set_LS%csf_set_LS%nwshells) then
         if(all_occupation(isubc).eq.0) then
            if(isubc.gt.1) then
               Li(isubc)  = 0;          L_i(isubc) = L_i(isubc-1)
               Si(isubc)  = 0;          S_i(isubc) = S_i(isubc-1)
               if(isubc .lt. asf_set_LS%csf_set_LS%nwshells) then
                  call setLS_job_count(isubc + 1, rez)
               else
                  if(ittk(S_i(isubc),L_i(isubc),J).eq.1)          &
                     call setLS_action(action_type, rez) !rez=rez+1
               end if
            else
               Li(isubc)  = 0;          L_i(isubc) = 0
               Si(isubc)  = 0;          S_i(isubc) = 0
            end if
         else
            N = all_occupation(isubc)
            l_shell = asf_set_LS%csf_set_LS%shell(isubc)%l
            call gettermLS(l_shell,N,LS_terms,number)
            do iterm=1, number, 1
!GG 2015.12.02 NIST
               N = all_occupation(isubc)
               l_shell = asf_set_LS%csf_set_LS%shell(isubc)%l
               call gettermLS(l_shell,N,LS_terms,numberGG)
!GG 2015.12.02 NIST
               Li(isubc)  = LS_terms(iterm)%LL
               Si(isubc)  = LS_terms(iterm)%S
               w(isubc)   = LS_terms(iterm)%w
               Q(isubc)   = LS_terms(iterm)%Q
               if(isubc.eq.1) then
                  L_i(isubc) = LS_terms(iterm)%LL
                  S_i(isubc) = LS_terms(iterm)%S
                  if(asf_set_LS%csf_set_LS%nwshells.gt.1) then
                     call setLS_job_count(isubc + 1, rez)
                  else
                     if(ittk(S_i(isubc),L_i(isubc),J).eq.1)          &
                        call setLS_action(action_type, rez) !rez=rez+1
                  end if
               else
                  tempLmax=L_i(isubc-1)+Li(isubc)
                  tempLmin=abs(L_i(isubc-1)-Li(isubc))
                  tempSmax=S_i(isubc-1)+Si(isubc)
                  tempSmin=abs(S_i(isubc-1)-Si(isubc))
                  do tempL = tempLmin, tempLmax, 2
                     L_i(isubc) = tempL
                     do tempS = tempSmin, tempSmax, 2
                        S_i(isubc) = tempS
                        if(isubc.lt.asf_set_LS%csf_set_LS%nwshells) then
                           call setLS_job_count(isubc+1,rez)
                        else
                           if(ittk(S_i(isubc),L_i(isubc),J).eq.1)        &
                            call setLS_action(action_type, rez) ! rez=rez+1
                        end if
                     end do   ! tempS
                  end do      ! tempL
               end if
            end do            ! iterm
         end if
      end if
      END SUBROUTINE setLS_job_count
!
!***********************************************************************
!                                                                      *
      FUNCTION setLS_equivalent_csfs(ncsf1,ncsf2)         result(rez)
!                                                                      *
!     This subroutine defines the "equivalency" of two csfs_jj         *
!     in the sence of generation of csfs_LS                            *
!     number of csfs_LS and corresponding quantum numbers.             *
!                                                                      *
!     Calls: itjpo, ispar, iq.                                         *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE ORB_C,           ONLY: NW
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE ispar_I
      USE IQ_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ncsf1, ncsf2
      logical rez
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isubc_LS, isubc_jj, NLS1, NLS2
!-----------------------------------------------
      rez= .true.
      IF(ITJPO(ncsf1) /= ITJPO(ncsf2)) THEN
         rez = .false.
      ELSE IF (ISPAR(ncsf1) /= ISPAR(ncsf2)) THEN
         rez = .false.
      ELSE
         do isubc_LS = asf_set_LS%csf_set_LS%nwcore + 1,  &
                       asf_set_LS%csf_set_LS%nwshells, 1
            NLS1 = 0;   NLS2 = 0
            do isubc_jj = 1, NW, 1
               if(nlLSval(isubc_LS) == nljjval(isubc_jj)) then
                   NLS1 = NLS1 + IQ(isubc_jj,ncsf1)
                   NLS2 = NLS2 + IQ(isubc_jj,ncsf2)
               end if
            end do
            if(NLS1.ne.NLS2) then
               rez=.false.
               return
            end if
         end do
      END IF
      END FUNCTION setLS_equivalent_csfs
!
      END SUBROUTINE setLS
!
!***********************************************************************
!                                                                      *
      FUNCTION traLSjj(jj_number,LS_number)          result(wa)
!                                                                      *
!     This procedure return the value of the transformation matrix     *
!     from jj- to LS-coupling scheme in the case of any number of      *
!     open shells.                                                     *
!                                                                      *
!     Calls: coefLSjj, iq, jqs, traLSjjmp.                             *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE ORB_C,           ONLY: NAK
      USE CONS_C,          ONLY: ZERO
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE IQ_I
      USE JQS_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jj_number,LS_number
      real(DOUBLE)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: total_number, shell_number, number_minus, number_plus, &
                 jj_minus, N_minus, Q_minus, J_1_minus,                 &
                 jj_plus, N_plus, Q_plus, J_1_plus, J_1,                &
                 l_shell, N_LS, W_1, Q_1, L_1, S_1
!-----------------------------------------------
      shell_number=asf_set_LS%csf_set_LS%nwcore+1
      number_minus=asf_set_LS%csf_set_LS%parent(shell_number)%parent_minus
      number_plus =asf_set_LS%csf_set_LS%parent(shell_number)%parent_plus
      if (number_minus+1 /= number_plus .and.  &
          number_minus*number_plus  /= 0) then
        stop "traLSjj(): program stop A."
      end if
!
      if (number_minus == 0) then
         jj_minus = iabs(NAK(number_plus))*2 - 3
         N_minus  = 0;    J_1_minus = 0;   Q_minus = (jj_minus + 1)/2
      else
         jj_minus  = iabs(NAK(number_minus))*2 - 1
         N_minus   = IQ(number_minus,jj_number)
         Q_minus   = (jj_minus +1)/2 - JQS(1,number_minus,jj_number)
         J_1_minus = JQS(3,number_minus,jj_number)-1
         J_1       = Jcoup(number_minus,jj_number)
      end if
!
      if (number_plus == 0) then
         jj_plus = iabs(NAK(number_minus))*2 + 1
         N_plus  = 0;    J_1_plus = 0;   Q_plus = (jj_plus + 1)/2
      else
         jj_plus  = iabs(NAK(number_plus))*2 - 1
         N_plus   = IQ(number_plus,jj_number)
         Q_plus   = (jj_plus + 1)/2 - JQS(1,number_plus,jj_number)
         J_1_plus = JQS(3,number_plus,jj_number)-1
         J_1       = Jcoup(number_plus,jj_number)
      end if
!
      l_shell=asf_set_LS%csf_set_LS%shell(shell_number)%l
      N_LS   =asf_set_LS%csf_set_LS%csf(LS_number)%occupation(shell_number)
      W_1    =asf_set_LS%csf_set_LS%csf(LS_number)%w(shell_number)
      Q_1    =2*l_shell+1-                                              &
              asf_set_LS%csf_set_LS%csf(LS_number)%seniority(shell_number)
!
      L_1     = asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number)
      S_1     = asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number)
!
      wa = coefLSjj(l_shell,N_LS,W_1,Q_1,L_1,S_1,J_1,                   &
                                  jj_minus,N_minus,Q_minus,J_1_minus,   &
                                  jj_plus,Q_plus,J_1_plus)
      if (abs(wa) > EPSNEW) then
         total_number =                                                 &
            asf_set_LS%csf_set_LS%nwshells - asf_set_LS%csf_set_LS%nwcore
         if (total_number == 1) then
            wa = wa
         else if (total_number >= 2) then
            do shell_number = asf_set_LS%csf_set_LS%nwcore + 2,         &
                              asf_set_LS%csf_set_LS%nwshells
               if (abs(wa) > EPSNEW) then
                  wa = wa * traLSjjmp(shell_number,jj_number,LS_number)
               end if
            end do
         else
            wa = zero
         end if
      else
         wa = zero
      end if
      END FUNCTION traLSjj
!
!***********************************************************************
!                                                                      *
      FUNCTION traLSjjmp(shell_number,jj_number,LS_number)   result(wa)
!                                                                      *
!     Return the value of main part of the transformation matrix from  *
!     jj- to LS-coupling scheme in the case of any number of open      *
!     shells.                                                          *
!                                                                      *
!     Calls: coefLSjj, iq, ixjtik, jqs, sixj, nine.                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE ORB_C,           ONLY: NAK
      USE CONS_C,          ONLY: ZERO, ONE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I
      USE IQ_I
      USE JQS_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: shell_number,jj_number,LS_number
      real(DOUBLE)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(DOUBLE) :: wa_sum, RAC6, RAC9
      integer      :: delta_J, number_minus, number_plus,   &
                      number_plus_1, number_minus_1,        &
                      jj_minus, N_minus, Q_minus,           &
                      jj_plus, N_plus, Q_plus,              &
                      J_i_min, J_i_max, J_i_minus,          &
                      J_i_plus, J_i, J_1_i, Jp_1_i, J_1_i1, &
                      l_shell, N_LS, W_i, Q_i, L_i, S_i,    &
                      L_1_i, S_1_i, L_1_i1, S_1_i1
!-----------------------------------------------
      wa     = zero;      wa_sum = zero
!
      number_minus = asf_set_LS%csf_set_LS%parent(shell_number)%parent_minus
      number_plus  = asf_set_LS%csf_set_LS%parent(shell_number)%parent_plus
      if (number_minus+1 /= number_plus .and.  &
          number_minus*number_plus  /= 0) then
         stop "tranLSjjmp(): program stop A."
      end if
!
      number_plus_1=asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_plus
      if (number_plus_1 == 0) then
         number_plus_1 =                                            &
                 asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_minus
      end if
!
      if (number_minus == 0) then
         jj_minus  = iabs(NAK(number_plus))*2 - 1
         N_minus   = 0;    J_i_minus = 0;    Q_minus = (jj_minus + 1)/2
         Jp_1_i    = Jcoup(number_plus_1,jj_number)
      else
         jj_minus  = iabs(NAK(number_minus))*2 - 1
         N_minus   = IQ(number_minus,jj_number)
         Q_minus   = (jj_minus + 1)/2 - JQS(1,number_minus,jj_number)
         J_i_minus = JQS(3,number_minus,jj_number)-1
         Jp_1_i    = Jcoup(number_minus,jj_number)
      end if
!
      if (number_plus == 0) then
         jj_plus  = iabs(NAK(number_minus))*2 - 1
         Q_plus   = (jj_plus + 1)/2
         N_plus   = 0;         J_i_plus = 0;         J_1_i   = Jp_1_i
      else
         jj_plus  = iabs(NAK(number_plus))*2 - 1
         N_plus   = IQ(number_plus,jj_number)
         Q_plus   = (jj_plus + 1)/2 - JQS(1,number_plus,jj_number)
         J_i_plus = JQS(3,number_plus,jj_number)-1
         J_1_i     = Jcoup(number_plus,jj_number)
      end if
!
      J_i_min = iabs(J_i_minus - J_i_plus);   J_i_max = J_i_minus + J_i_plus
      !
      l_shell = asf_set_LS%csf_set_LS%shell(shell_number)%l
      N_LS    = asf_set_LS%csf_set_LS%csf(LS_number)%occupation(shell_number)
      W_i     = asf_set_LS%csf_set_LS%csf(LS_number)%w(shell_number)
      Q_i     = 2 * l_shell + 1 -                                          &
                asf_set_LS%csf_set_LS%csf(LS_number)%seniority(shell_number)
      if (N_LS == N_minus + N_plus) then
         L_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number)
         S_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number)
         !
         L_1_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellLX(shell_number)
         S_1_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellSX(shell_number)
         if (shell_number == 2) then
            number_minus_1  =                                           &
               asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_minus
            if (number_minus_1 == 0) then
               number_minus_1 =                                         &
               asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_plus
            end if
            J_1_i1 =JQS(3,number_minus_1,jj_number)-1
            L_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number-1)
            S_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number-1)
         else
            J_1_i1 =Jcoup(number_plus_1,jj_number)
            L_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellLX(shell_number-1)
            S_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellSX(shell_number-1)
         end if
!
         do J_i = J_i_min, J_i_max, 2
            delta_J = ixjtik(J_i_minus,J_i_plus,J_i,J_1_i,J_1_i1,Jp_1_i)
            if (delta_J /= 0) then
               call nine(L_1_i1,S_1_i1,J_1_i1,L_i,S_i,J_i,L_1_i,S_1_i,J_1_i,  &
                      1,delta_J,RAC9)
               if (delta_J /= 0) then
                 call sixj(J_i_minus,J_i_plus,J_i,J_1_i,J_1_i1,Jp_1_i,0,RAC6)
                 call nine(L_1_i1,S_1_i1,J_1_i1,L_i,S_i,J_i,L_1_i,S_1_i,J_1_i,&
                      0,delta_J,RAC9)
                 wa_sum = wa_sum + (J_i + one) * RAC6 * RAC9 *                &
                      coefLSjj(l_shell,N_LS,W_i,Q_i,L_i,S_i,J_i,              &
                                       jj_minus,N_minus,Q_minus,J_i_minus,    &
                                       jj_plus,Q_plus,J_i_plus)
               end if
            end if
         end do
         wa = wa_sum * sqrt((Jp_1_i+one)*(J_1_i1+one)*(L_1_i+one)*(S_1_i+one))
         if (mod(J_i_minus+J_i_plus+J_1_i1+J_1_i,4) /= 0) wa = - wa
      else
         wa = zero
      end if
      END FUNCTION traLSjjmp
!
!***********************************************************************
!                                                                      *
      SUBROUTINE uniquelsj
!                                                                      *
!     Subroutine defines a unique labels for energy levels             *
!                                                                      *
!     Written by  G. Gaigalas                        Vilnius, May 2017 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     The maximum number of levels in the list
      INTEGER, PARAMETER :: Lev_No  = 100
!     The maximum number of mixing coefficients of ASF expension
      INTEGER, PARAMETER :: Vec_No  = 100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER, DIMENSION(Vec_No,Lev_No)    :: COUPLING*256
      REAL(DOUBLE), DIMENSION(Vec_No,Lev_No) :: COMP, MIX
      CHARACTER, DIMENSION(Lev_No)           :: Str_No*3,Str_J*5,Str_P*1
      CHARACTER, DIMENSION(Lev_No)           :: OPT_COUPLING*256
      CHARACTER, DIMENSION(Vec_No)           :: tmp_COUPLING*256
      REAL(DOUBLE), DIMENSION(Lev_No)        :: ENERGY, PRO
      REAL(DOUBLE), DIMENSION(Vec_No)        :: tmp_COMP, tmp_MIX
      INTEGER, DIMENSION(Lev_No)             :: MAS_MAX, ICOUNT, IPRGG
      CHARACTER    :: RECORD*53
      CHARACTER    :: FORM*11
      CHARACTER    :: Str*7
      REAL(DOUBLE) :: MAX_COMP
      INTEGER      :: IOS, IERR, INUM, NUM_Lev, I1,I2,I3,I4
      INTEGER      :: IPR, IPRF, NUM_OPT, Lev_OPT, New_J
      INTEGER*4    :: last, length
!-----------------------------------------------
      REWIND(57)
      read (57, '(1A53)', IOSTAT=IOS) RECORD
      write(73, '(1A53)') RECORD
!
      DO WHILE (New_J < 1)
        NUM_Lev = 0
        New_J = 1
        DO
          if(NUM_Lev > Lev_No) then
            print*, "Please extand the arrays. Now we have Lev_No=",   &
                                                                  Lev_No
            stop
          end if
          NUM_Lev = NUM_Lev + 1
          read(57,'(A3,A5,5X,A1,8X,F16.9,5X,F7.3)')Str_No(NUM_Lev),    &
              Str_J(NUM_Lev),Str_P(NUM_Lev),Energy(NUM_Lev),PRO(NUM_Lev)
          INUM = 0
          DO
            INUM = INUM + 1
            if(INUM > Vec_No) then
              print*, "Please extand the arrays. Now we have Vec_No=", &
                                                                  Vec_No
              stop
              end if
            read(57,'(A7,F12.8,3X,F11.8,3X,A)',IOSTAT=IOS)Str,         &
             MIX(INUM,NUM_Lev),COMP(INUM,NUM_Lev),COUPLING(INUM,NUM_Lev)
            if (IOS==-1) then
               MAS_MAX(NUM_Lev) = INUM - 1
               go to 1
            endif
            if(Str /= '       ') then
               backspace(57)
               MAS_MAX(NUM_Lev) = INUM - 1
             exit
            else if                                                    &
              (MIX(INUM,NUM_Lev)==0.00.and.COMP(INUM,NUM_Lev)==0.00)then
               MAS_MAX(NUM_Lev) = INUM - 1
               New_J = 0
               go to 1
            end if
          END DO
        END DO
   1    CONTINUE
!
        NUM_OPT = 0
        ICOUNT = 1
        IPRGG = 1
        write(72,'(A)')                                                &
                    "          Composition  Serial No.         Coupling"
        write(72,'(A)') "                       of compos."
        write(72,'(A5,A)') " J = ",trim(Str_J(NUM_Lev))
        write(72,'(A)')                                                &
                    "--------------------------------------------------"
        DO I1 = 1, NUM_Lev
   2      MAX_COMP = 0.0
          DO I2 = 1, NUM_Lev
            if(ICOUNT(I2) == 0) CYCLE
            if(COMP(1,I2) > MAX_COMP) then
               Lev_OPT = I2
               MAX_COMP = COMP(1,I2)
            end if
          END DO
          NUM_OPT = NUM_OPT + 1
          ICOUNT(Lev_OPT) = 0
          IPR = 1
          if(NUM_OPT == 1) then
            IPR = 1
          else
            I3 = 1
            DO WHILE (I3 < NUM_OPT)
             if(trim(OPT_COUPLING(I3))==trim(COUPLING(IPR,Lev_OPT)))then
               IPR = IPR + 1
               if(IPR > Vec_No) then
                 print*,"Please extand the arrays. Now we have Vec_No="&
                                                                 ,Vec_No
                 stop
               end if
               I3 = 1
             else
               I3 = I3 + 1
             end if
            END DO
          end if
          IPRGG(Lev_OPT) = IPR + IPRGG(Lev_OPT) - 1
          if(IPRGG(Lev_OPT) > MAS_MAX(Lev_OPT)                        &
                                        .AND. IPRGG(Lev_OPT) /= 1) then
             print*,                                                   &
             "The program is not able perform the identification for", &
             " level = ",Lev_OPT
             stop
          end if
          if(IPR == 1) then
            OPT_COUPLING(NUM_OPT) = COUPLING(IPR,Lev_OPT)
            IPRF = IPRGG(Lev_OPT)
            write(72,'(A,I4,2X,F12.9,I5,3X,A,A,I4)')"Pos",Lev_OPT,     &
            COMP(IPR,Lev_OPT),IPRF,trim(OPT_COUPLING(NUM_OPT))
            IF(IPRF > 1) THEN
              tmp_MIX = 0
              tmp_COMP = 0
              tmp_COUPLING = ""
              I4 = MAS_MAX(Lev_OPT)
              tmp_MIX(2:IPRF)          = MIX(I4-IPRF+2:I4,Lev_OPT)
              tmp_MIX(IPRF+1:I4)       = MIX(2:I4-IPRF+1,Lev_OPT)
              MIX(2:I4,Lev_OPT)        = tmp_MIX(2:I4)
!
              tmp_COMP(2:IPRF)         = COMP(I4-IPRF+2:I4,Lev_OPT)
              tmp_COMP(IPRF+1:I4)      = COMP(2:I4-IPRF+1,Lev_OPT)
              COMP(2:I4,Lev_OPT)       = tmp_COMP(2:I4)
!
              tmp_COUPLING(2:IPRF)    = COUPLING(I4-IPRF+2:I4,Lev_OPT)
              tmp_COUPLING(IPRF+1:I4) = COUPLING(2:I4-IPRF+1,Lev_OPT)
              COUPLING(2:I4,Lev_OPT)  = tmp_COUPLING(2:I4)
            END IF
          end if
          if(IPR > 1) then
            tmp_MIX = 0
            tmp_COMP = 0
            tmp_COUPLING = ""
            I4 = MAS_MAX(Lev_OPT)
            tmp_MIX(1:1) = MIX(IPR:IPR,Lev_OPT)
            tmp_MIX(I4-IPR+2:I4)      = MIX(1:IPR-1,Lev_OPT)
            tmp_MIX(2:I4-IPR+1)       = MIX(IPR+1:I4,Lev_OPT)
            MIX(1:I4,Lev_OPT)         = tmp_MIX(1:I4)
!
            tmp_COMP(1:1)             = COMP(IPR:IPR,Lev_OPT)
            tmp_COMP(I4-IPR+2:I4)     = COMP(1:IPR-1,Lev_OPT)
            tmp_COMP(2:I4-IPR+1)      = COMP(IPR+1:I4,Lev_OPT)
            COMP(1:I4,Lev_OPT)        = tmp_COMP(1:I4)
!
            tmp_COUPLING(1:1)         = COUPLING(IPR:IPR,Lev_OPT)
            tmp_COUPLING(I4-IPR+2:I4) = COUPLING(1:IPR-1,Lev_OPT)
            tmp_COUPLING(2:I4-IPR+1)  = COUPLING(IPR+1:I4,Lev_OPT)
            COUPLING(1:I4,Lev_OPT)    = tmp_COUPLING(1:I4)
!
            ICOUNT(Lev_OPT) = 1
            NUM_OPT = NUM_OPT - 1
            go to 2
          end if
        END DO
        write(72,'(A)')                                                &
                    "--------------------------------------------------"
        write(72,'(A)')""
!
        DO I1 = 1, NUM_Lev
          write(73,'(A3,A5,5X,A1,8X,F16.9,5X,F7.3,A)')                 &
          Str_No(I1),Str_J(I1),Str_P(I1),Energy(I1),PRO(I1),"%"
          DO I2 = 1, MAS_MAX(I1)
            write(73,'(A7,F12.8,3X,F11.8,3X,A)')                       &
            Str,MIX(I2,I1),COMP(I2,I1),trim(COUPLING(I2,I1))
          END DO
        END DO
        if(New_J < 1) write(73,'(A)') " "
!
      end do
      END SUBROUTINE uniquelsj
!
END MODULE jj2lsj_code
