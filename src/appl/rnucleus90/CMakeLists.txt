add_executable(rnucleus
    nucleus_m.f90
    skfun.f90
    skfun_I.f90
    estrms.f90
    estrms_I.f90
    getcpr.f90
    getcpr_I.f90
    geniso.f90
)
target_link_libraries_Fortran(rnucleus PUBLIC 9290 mod)
install(TARGETS rnucleus DESTINATION bin/)
