add_executable(rhfszeeman95
    opt6_C.f90
    engouth_I.f90
    gethfd_I.f90
    hfszeeman_I.f90
    matelt_I.f90
    rinthf_I.f90
    setdbg_I.f90
    setsum_I.f90
    strsum_I.f90
    getmixblock_I.f90
    engouth.f90
    gethfd.f90
    hfszeeman.f90
    hfszeeman06.f90
    matelt.f90
    rinthf.f90
    setdbg.f90
    setsum.f90
    strsum.f90
    getmixblock.f90
)
target_link_libraries_Fortran(rhfszeeman95 PUBLIC 9290 mod rang90 mcp90 9290 mod)
install(TARGETS rhfszeeman95 DESTINATION bin/)
