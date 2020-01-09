add_executable(rangular_mpi
    fndbeg.f90
    fndbeg_I.f90
    getinf.f90
    getinf_I.f90
    outsdampi.f90
    outsdampi_I.f90
    setdbg.f90
    setdbg_I.f90
    setdbgmpi.f90
    setdbgmpi_I.f90
    setmcp.f90
    setmcp_I.f90
    setmcpmpi.f90
    setmcpmpi_I.f90
    setsda.f90
    setsda_I.f90
    setsum.f90
    setsum_I.f90
    settmp.f90
    settmp_I.f90
    sort.f90
    sort_I.f90
    strsum.f90
    strsum_I.f90
    mcpmpi_gg.f90
    mcpmpi_gg_I.f90
    genmcpmpi.f90
)
target_link_libraries_Fortran(rangular_mpi PUBLIC rang90 mcp90 mpiu90 9290 mod)
target_link_libraries(rangular_mpi PRIVATE ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
target_link_libraries(rangular_mpi PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
target_include_directories(rangular_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(rangular_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
set_target_properties(rangular_mpi PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
install(TARGETS rangular_mpi DESTINATION bin/)
