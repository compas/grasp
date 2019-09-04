add_executable(rnucleus
    estrms.f90
    estrms_I.f90
    geniso.f90
    getcpr.f90
    getcpr_I.f90
    skfun.f90
    skfun_I.f90
)
target_link_libraries_Fortran(rnucleus PUBLIC mod 9290)
install(TARGETS rnucleus DESTINATION bin/)
