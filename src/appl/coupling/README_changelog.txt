The following changes was made after publication in CPC:
----------------------------------------------
module Coupling_main.f90
       subroutine print_DataBasis

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_evaluation.f90
       SUBROUTINE getchLS
       SUBROUTINE getchjj

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_inside_shell.f90
       FUNCTION gettermLSQ
       subroutine form_csfs_cLSJ3(itype)
       subroutine  define_number_of_csfs_cLSJ3(irez)
       subroutine matrix_LS_cLSJ3(icsf_LS,icsf_cLSJ3,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LScLSJ3.f90
       subroutine form_csfs_cLSJ3(itype)
       subroutine  define_number_of_csfs_cLSJ3(irez)
       subroutine matrix_LS_cLSJ3(icsf_LS,icsf_cLSJ3,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LSjj1.f90
       subroutine form_csfs_jj1(itype)
       subroutine  define_number_of_csfs_jj1(irez)
       subroutine matrix_LS_jj1(icsf_LS,icsf_jj1,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LSjj2.f90
       subroutine form_csfs_jj2(itype)
       subroutine  define_number_of_csfs_jj2(irez)
       subroutine matrix_LS_jj2(icsf_LS,icsf_jj2,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LSjj2.f90
       subroutine form_csfs_jj2(itype)
       subroutine  define_number_of_csfs_jj2(irez)
       subroutine matrix_LS_jj2(icsf_LS,icsf_jj2,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LSjj3.f90
       subroutine form_csfs_jj3(itype)
       subroutine  define_number_of_csfs_jj3(irez)
       subroutine matrix_LS_jj3(icsf_LS,icsf_jj3,rez)

 G.G.       2020.02.04.
----------------------------------------------
module Coupling_transform_cfg_LSLScjj.f90
       subroutine form_csfs_LScjj(itype)
       subroutine  define_number_of_csfs_LScjj(irez)
       subroutine matrix_LS_LScjj(icsf_LS,icsf_LScjj,rez)

 G.G.       2020.02.04.
----------------------------------------------
module sqlsf.f
       CHARACTER*4 FUNCTION JVAL(IVALUE)

 G.G.       2020.02.05.
----------------------------------------------
module Coupling.f90
       subroutine open_files
       subroutine get_parameters

 G.G.       2020.02.05.
----------------------------------------------
module Coupling_transform_cfg_LSJJ.f90
       subroutine form_csfs_JJ(itype)

 G.G.       2020.02.11.
----------------------------------------------
module Coupling_transform_cfg_JJJK.f90
       subroutine form_csfs_JK(itype)

 G.G.       2020.02.11.
----------------------------------------------
module Coupling_transform_cfg_LSLK.f90
       subroutine form_csfs_LK(itype)

 G.G.       2020.02.11.
----------------------------------------------
module Coupling_transform_cfg_LSLK3.f90
       subroutine form_csfs_LK3(itype)

 G.G.       2020.02.11.
----------------------------------------------
module Coupling_transform_cfg_LSLS3.f90
       subroutine form_csfs_LS3(itype)

 G.G.       2020.02.11.
----------------------------------------------
