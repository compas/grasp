add_library(mcp90 STATIC
    cxk.f90
    cxk_I.f90
    talk.f90
    talk_I.f90
)
setup_fortran_modules(mcp90)
target_link_libraries_Fortran(mcp90 PRIVATE mod 9290)
install(TARGETS mcp90 DESTINATION lib/)
