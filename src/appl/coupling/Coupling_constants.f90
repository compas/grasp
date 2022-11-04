!
!***********************************************************************
!                                                                      *
      module Coupling_constants
!                                                                      *
!     Written by G. Gaigalas,                                          *
!                                                                      *
!     NIST                                     last update: Oct 2015   *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter, public :: i4b = selected_int_kind(9)
      integer, parameter, public :: i2b = selected_int_kind(4)
      integer, parameter, public :: i1b = selected_int_kind(2)
      integer, parameter, public :: sp  = kind(1.0)
      integer, parameter, public :: dp  =                              &
                                 selected_real_kind(2*precision(1.0_sp))
      real(kind=dp), parameter, public :: ZERO_dp = 0.0
      real(kind=dp), parameter, public :: ONE_dp = 1.0
      real(kind=dp), parameter, public :: TWO_dp = 2.0
      real(kind=dp), public ::    dp_coeff_precision     = 0.00000001
      real(kind=dp), public ::    dp_coeff_precision_t10 = 0.000001
      real(kind=sp), public ::    sp_coeff_precision     = 0.0000001
!     read unit
      integer, parameter, public :: iread_csf   = 1        !for csf (usually file "csf.inp")
      integer, parameter, public :: iread_coefs = 2        !for J-matrices and coefficients (usially file "name.j")
      integer, parameter, public :: iwrite_log  = 3
      integer, parameter, public :: iwrite_log2 = 3
      integer, parameter, public :: iwrite_expansions           = 16
      integer, parameter, public :: iwrite_cfg_expansions_JJ    =  7
      integer, parameter, public :: iwrite_cfg_expansions_LK    =  8
      integer, parameter, public :: iwrite_cfg_expansions_JK    =  9
      integer, parameter, public :: iwrite_cfg_expansions_LS3   = 12
      integer, parameter, public :: iwrite_cfg_expansions_LSJ3  = 13
      integer, parameter, public :: iwrite_cfg_expansions_LK3   = 14
      integer, parameter, public :: iwrite_cfg_expansions_JK3   = 15
      integer, parameter, public :: iwrite_cfg_expansions_cLSJ3 = 17
      integer, parameter, public :: iwrite_cfg_expansions_LScjj = 18
      integer, parameter, public :: iwrite_cfg_expansions_jj123 = 19
      integer, parameter, public :: iwrite_classifications      = 10
      integer, parameter, public :: iwrite_DataBase_LS          = 50
      integer, parameter, public :: iwrite_DataBase_JJ          = 51
      integer, parameter, public :: iwrite_DataBase_LK          = 52
      integer, parameter, public :: iwrite_DataBase_JK          = 53
      integer, parameter, public :: iwrite_DataBase_LS3         = 64
      integer, parameter, public :: iwrite_DataBase_LSJ3        = 65
      integer, parameter, public :: iwrite_DataBase_LK3         = 56
      integer, parameter, public :: iwrite_DataBase_JK3         = 57
      integer, parameter, public :: iwrite_DataBase_cLSJ3       = 58
      integer, parameter, public :: iwrite_DataBase_LScjj       = 59
      integer, parameter, public :: iwrite_DataBase_jj123       = 60
      integer, parameter, public ::                              &
                          iwrite_classifications_evaluations =  11
      integer, parameter, public :: iwrite_suggestions       =  70
!
      integer, parameter, public :: NOTATION =1  ! for ASD clasification
!      integer, parameter, public :: NOTATION =0  ! for ATSP2K and GRASP2K
!-----------------------------------------------
      end module Coupling_constants
