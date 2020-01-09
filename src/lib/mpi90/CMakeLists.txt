add_library(mpiu90 STATIC
    mpi_C.f90
    mpiu.f90
    cpath.f90
    cpath_I.f90
    cslhmpi.f90
    cslhmpi_I.f90
    iniestmpi.f90
    iniestmpi_I.f90
    lodcslmpi.f90
    lodcslmpi_I.f90
    lodrwfmpi.f90
    lodrwfmpi_I.f90
    setisompi.f90
    setisompi_I.f90
    setrwfmpi.f90
    setrwfmpi_I.f90
    spicmvmpi.f90
    spicmvmpi_I.f90
    sys_chdir.f90
    sys_chdir_I.f90
    sys_getwd.f90
    sys_getwd_I.f90
    sys_mkdir.f90
    sys_mkdir_I.f90
)
setup_fortran_modules(mpiu90)
target_link_libraries_Fortran(mpiu90 PRIVATE mod 9290)
target_include_directories(mpiu90 PRIVATE ${MPI_Fortran_INCLUDE_PATH})"
target_link_libraries(mpiu90 PRIVATE ${MPI_Fortran_LIBRARIES})"
set_target_properties(mpiu90 PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
