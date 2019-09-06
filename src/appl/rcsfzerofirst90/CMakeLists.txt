add_executable(rcsfzerofirst
    RCSFzerofirst.f90

    lodcsl_Part.f90 lodcsl_Part_I.f90
    lodcsl_Zero.f90 lodcsl_Zero_I.f90
    set_CSF_ZFlist.f90 set_CSF_ZFlist_I.f90
    set_CSF_number.f90 set_CSF_number_I.f90
)
target_link_libraries_Fortran(rcsfzerofirst PUBLIC mod 9290)
install(TARGETS rcsfzerofirst DESTINATION bin/)
